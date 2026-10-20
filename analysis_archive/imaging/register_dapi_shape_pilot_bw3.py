from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy.optimize import minimize
from skimage import exposure, filters, measure, morphology, transform


ROOT = Path(r"data/external/local_77fbce6c2be7")
OME = ROOT / "01_ome_tiff"
OUTPUT = ROOT / "02_registration_and_lanes" / "00_audit"
OUTPUT.mkdir(parents=True, exist_ok=True)
FIXED_PATH = OME / "KC_Mig_浅ⅡBW3.ome.tiff"
MOVING_PATH = OME / "Fib_K14_BW3 PDGFRα 520+P63 570+KLF4 620+KRT14 690(1).ome.tiff"
LEVEL = 7


def load_dapi(path: Path) -> np.ndarray:
    with tifffile.TiffFile(path) as tif:
        level = tif.series[0].levels[min(LEVEL, len(tif.series[0].levels) - 1)]
        return np.asarray(level.asarray()[0], dtype=np.uint8)


def normalize(array: np.ndarray) -> np.ndarray:
    positive = array[array > 0]
    if positive.size == 0:
        return np.zeros(array.shape, dtype=np.float32)
    lo, hi = np.percentile(positive, (0.5, 99.8))
    if hi <= lo:
        hi = lo + 1
    scaled = np.clip((array.astype(np.float32) - lo) / (hi - lo), 0, 1)
    return exposure.equalize_adapthist(scaled, kernel_size=32, clip_limit=0.02).astype(np.float32)


def tissue_mask(image: np.ndarray) -> np.ndarray:
    smooth = filters.gaussian(image, sigma=3.0, preserve_range=True)
    positive = smooth[smooth > 0.01]
    threshold = filters.threshold_otsu(positive) if positive.size else 1.0
    mask = smooth > max(0.035, threshold * 0.55)
    mask = morphology.binary_closing(mask, morphology.disk(5))
    mask = morphology.binary_dilation(mask, morphology.disk(3))
    mask = morphology.remove_small_objects(mask, min_size=200)
    mask = morphology.remove_small_holes(mask, area_threshold=1000)
    labels = measure.label(mask)
    regions = sorted(measure.regionprops(labels), key=lambda region: region.area, reverse=True)
    if regions:
        keep = {region.label for region in regions if region.area >= regions[0].area * 0.02}
        mask = np.isin(labels, list(keep))
    return mask


def orientation_matrix(shape: tuple[int, int], rotation: int, flip_lr: bool) -> tuple[np.ndarray, tuple[int, int]]:
    h, w = shape
    if rotation == 0:
        matrix, out_shape = np.eye(3), (h, w)
    elif rotation == 90:
        matrix = np.array([[0, 1, 0], [-1, 0, w - 1], [0, 0, 1]], dtype=float)
        out_shape = (w, h)
    elif rotation == 180:
        matrix = np.array([[-1, 0, w - 1], [0, -1, h - 1], [0, 0, 1]], dtype=float)
        out_shape = (h, w)
    elif rotation == 270:
        matrix = np.array([[0, -1, h - 1], [1, 0, 0], [0, 0, 1]], dtype=float)
        out_shape = (w, h)
    else:
        raise ValueError(rotation)
    if flip_lr:
        matrix = np.array([[-1, 0, out_shape[1] - 1], [0, 1, 0], [0, 0, 1]], dtype=float) @ matrix
    return matrix, out_shape


def bbox_and_center(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    yy, xx = np.nonzero(mask)
    if xx.size == 0:
        raise RuntimeError("Empty tissue mask")
    bbox = np.array([xx.min(), yy.min(), xx.max() + 1, yy.max() + 1], dtype=float)
    center = np.array([np.median(xx), np.median(yy)], dtype=float)
    return bbox, center


def dice(a: np.ndarray, b: np.ndarray) -> float:
    denom = int(a.sum() + b.sum())
    return 2.0 * int(np.logical_and(a, b).sum()) / denom if denom else 0.0


fixed = normalize(load_dapi(FIXED_PATH))
moving = normalize(load_dapi(MOVING_PATH))
fixed_mask = tissue_mask(fixed)
moving_mask = tissue_mask(moving)
fixed_bbox, fixed_center = bbox_and_center(fixed_mask)
fixed_box_size = fixed_bbox[2:] - fixed_bbox[:2]
fixed_canvas_center = np.array([(fixed.shape[1] - 1) / 2, (fixed.shape[0] - 1) / 2])

reports: list[dict] = []
best: dict | None = None

for rotation in (0, 90, 180, 270):
    for flip_lr in (False, True):
        orientation, oriented_shape = orientation_matrix(moving.shape, rotation, flip_lr)
        oriented_mask = transform.warp(
            moving_mask.astype(np.uint8),
            inverse_map=transform.AffineTransform(matrix=orientation).inverse,
            output_shape=oriented_shape,
            order=0,
            preserve_range=True,
        ) > 0
        moving_bbox, moving_center = bbox_and_center(oriented_mask)
        moving_box_size = moving_bbox[2:] - moving_bbox[:2]
        sx, sy = np.clip(fixed_box_size / moving_box_size, 0.75, 1.35)
        initial = np.array(
            [[sx, 0, fixed_center[0] - sx * moving_center[0]], [0, sy, fixed_center[1] - sy * moving_center[1]], [0, 0, 1]],
            dtype=float,
        ) @ orientation

        def delta_matrix(parameters: np.ndarray) -> np.ndarray:
            angle, log_sx, log_sy, shear, tx, ty = parameters
            theta = np.deg2rad(angle)
            scale_x = np.exp(log_sx)
            scale_y = np.exp(log_sy)
            linear = np.array(
                [
                    [np.cos(theta) * scale_x, (-np.sin(theta) + shear * np.cos(theta)) * scale_y, 0],
                    [np.sin(theta) * scale_x, (np.cos(theta) + shear * np.sin(theta)) * scale_y, 0],
                    [0, 0, 1],
                ]
            )
            to_origin = np.array([[1, 0, -fixed_canvas_center[0]], [0, 1, -fixed_canvas_center[1]], [0, 0, 1]])
            back = np.array([[1, 0, fixed_canvas_center[0] + tx], [0, 1, fixed_canvas_center[1] + ty], [0, 0, 1]])
            return back @ linear @ to_origin

        def objective(parameters: np.ndarray) -> float:
            total = delta_matrix(parameters) @ initial
            model = transform.AffineTransform(matrix=total)
            warped_mask = transform.warp(
                moving_mask.astype(np.uint8), model.inverse, output_shape=fixed.shape, order=0, preserve_range=True
            ) > 0
            dice_value = dice(fixed_mask, warped_mask)
            overlap = np.logical_and(fixed_mask, warped_mask)
            if overlap.sum() < 100:
                return 10.0
            warped_image = transform.warp(moving, model.inverse, output_shape=fixed.shape, order=1, preserve_range=True)
            f = fixed[overlap]
            m = warped_image[overlap]
            ncc = float(np.corrcoef(f, m)[0, 1]) if f.std() > 0 and m.std() > 0 else 0.0
            return -(0.82 * dice_value + 0.18 * max(-1.0, ncc))

        bounds = [(-12, 12), (-0.18, 0.18), (-0.18, 0.18), (-0.16, 0.16), (-60, 60), (-60, 60)]
        optimized = minimize(
            objective,
            np.zeros(6),
            method="Powell",
            bounds=bounds,
            options={"maxiter": 180, "xtol": 0.03, "ftol": 0.0002},
        )
        total = delta_matrix(optimized.x) @ initial
        model = transform.AffineTransform(matrix=total)
        warped_mask = transform.warp(
            moving_mask.astype(np.uint8), model.inverse, output_shape=fixed.shape, order=0, preserve_range=True
        ) > 0
        warped = transform.warp(moving, model.inverse, output_shape=fixed.shape, order=1, preserve_range=True)
        overlap = np.logical_and(fixed_mask, warped_mask)
        dice_value = dice(fixed_mask, warped_mask)
        if overlap.sum() > 100:
            f = fixed[overlap]
            m = warped[overlap]
            ncc = float(np.corrcoef(f, m)[0, 1]) if f.std() > 0 and m.std() > 0 else 0.0
        else:
            ncc = 0.0
        score = 0.82 * dice_value + 0.18 * ncc
        report = {
            "rotation_degrees_ccw": rotation,
            "flip_left_right": flip_lr,
            "tissue_dice": dice_value,
            "overlap_ncc": ncc,
            "score": score,
            "optimizer_success": bool(optimized.success),
            "optimizer_parameters": optimized.x.tolist(),
            "matrix_moving_level7_to_fixed_level7": total.tolist(),
        }
        reports.append(report)
        if best is None or score > best["score"]:
            best = {**report, "warped": warped, "warped_mask": warped_mask}

assert best is not None
warped = best.pop("warped")
best.pop("warped_mask")
overlay = np.stack(
    [np.clip(fixed * 255, 0, 255), np.clip(warped * 255, 0, 255), np.zeros_like(fixed)], axis=-1
).astype(np.uint8)
Image.fromarray(overlay).save(OUTPUT / "DAPI_shape_registration_best.png")
Image.fromarray((fixed_mask.astype(np.uint8) * 255)).save(OUTPUT / "DAPI_fixed_tissue_mask.png")
(OUTPUT / "DAPI_shape_registration_report.json").write_text(
    json.dumps({"best": best, "candidates": reports}, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(best, ensure_ascii=False, indent=2))
