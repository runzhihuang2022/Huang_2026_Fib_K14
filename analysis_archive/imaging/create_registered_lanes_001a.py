from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image
from scipy.optimize import minimize
from skimage import exposure, filters, measure, morphology, transform


ROOT = Path(r"data/external/local_0eb0c38c2043")
INPUT_DIR = ROOT / "01_ome_tiff"
OUTPUT_DIR = ROOT / "02_registration_and_lanes"
OME = INPUT_DIR
OUTPUT = OUTPUT_DIR
AUDIT = OUTPUT / "00_audit"
ALIGNED = OUTPUT / "01_aligned_channels"
LANES = OUTPUT / "02_lanes"
for directory in (AUDIT, ALIGNED, LANES):
    directory.mkdir(parents=True, exist_ok=True)

FIXED_PATH = OME / "KC_Mig_编号001A.ome.tiff"
MOVING_PATH = OME / "Fib_K14_编号001A PDGFRα 520+P63 570+KLF4 620+KRT14 690.ome.tiff"
COARSE_LEVEL = 7
# Refine at the same audited resolution as the orientation pilot.
REFINE_LEVEL = 7
OUTPUT_LEVEL = 2
BASE_PIXEL_SIZE_UM = 0.25284451246261597


def load_channel(path: Path, level_index: int, channel_index: int) -> tuple[np.ndarray, tuple[int, int]]:
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        level_index = min(level_index, len(series.levels) - 1)
        level = series.levels[level_index]
        array = np.asarray(level.asarray(key=channel_index), dtype=np.uint8)
        base_shape = tuple(series.shape[-2:])
    return array, base_shape


def normalize_registration(array: np.ndarray) -> np.ndarray:
    positive = array[array > 0]
    if positive.size == 0:
        return np.zeros(array.shape, dtype=np.float32)
    lo, hi = np.percentile(positive, (0.5, 99.8))
    if hi <= lo:
        hi = lo + 1.0
    scaled = np.clip((array.astype(np.float32) - lo) / (hi - lo), 0, 1)
    return exposure.equalize_adapthist(scaled, kernel_size=48, clip_limit=0.018).astype(np.float32)


def normalize_display(array: np.ndarray, low_percentile: float = 1.0, high_percentile: float = 99.8) -> np.ndarray:
    positive = array[array > 0]
    if positive.size == 0:
        return np.zeros(array.shape, dtype=np.uint8)
    lo, hi = np.percentile(positive, (low_percentile, high_percentile))
    if hi <= lo:
        hi = lo + 1.0
    return np.clip((array.astype(np.float32) - lo) * 255.0 / (hi - lo), 0, 255).astype(np.uint8)


def tissue_mask(image: np.ndarray) -> np.ndarray:
    smooth = filters.gaussian(image, sigma=4.0, preserve_range=True)
    positive = smooth[smooth > 0.01]
    threshold = filters.threshold_otsu(positive) if positive.size else 1.0
    mask = smooth > max(0.03, threshold * 0.55)
    mask = morphology.binary_closing(mask, morphology.disk(7))
    mask = morphology.binary_dilation(mask, morphology.disk(4))
    mask = morphology.remove_small_objects(mask, min_size=500)
    mask = morphology.remove_small_holes(mask, area_threshold=3000)
    labels = measure.label(mask)
    regions = sorted(measure.regionprops(labels), key=lambda region: region.area, reverse=True)
    if regions:
        keep = {region.label for region in regions if region.area >= regions[0].area * 0.02}
        mask = np.isin(labels, list(keep))
    return mask


def scale_matrix_between_levels(
    matrix: np.ndarray,
    moving_from_shape: tuple[int, int],
    fixed_from_shape: tuple[int, int],
    moving_to_shape: tuple[int, int],
    fixed_to_shape: tuple[int, int],
) -> np.ndarray:
    moving_to_from = np.diag(
        [moving_from_shape[1] / moving_to_shape[1], moving_from_shape[0] / moving_to_shape[0], 1.0]
    )
    fixed_from_to = np.diag(
        [fixed_to_shape[1] / fixed_from_shape[1], fixed_to_shape[0] / fixed_from_shape[0], 1.0]
    )
    return fixed_from_to @ matrix @ moving_to_from


def dice(a: np.ndarray, b: np.ndarray) -> float:
    denominator = int(a.sum() + b.sum())
    return 2.0 * int(np.logical_and(a, b).sum()) / denominator if denominator else 0.0


print(f"Reading inputs from: {INPUT_DIR}", flush=True)
coarse_report = json.loads((AUDIT / "DAPI_shape_registration_report.json").read_text(encoding="utf-8"))["best"]
coarse_matrix = np.asarray(coarse_report["matrix_moving_level7_to_fixed_level7"], dtype=float)

fixed_refine_raw, _ = load_channel(FIXED_PATH, REFINE_LEVEL, 0)
moving_refine_raw, _ = load_channel(MOVING_PATH, REFINE_LEVEL, 0)
fixed_coarse_raw, _ = load_channel(FIXED_PATH, COARSE_LEVEL, 0)
moving_coarse_raw, _ = load_channel(MOVING_PATH, COARSE_LEVEL, 0)
initial_refine = scale_matrix_between_levels(
    coarse_matrix,
    moving_coarse_raw.shape,
    fixed_coarse_raw.shape,
    moving_refine_raw.shape,
    fixed_refine_raw.shape,
)
fixed_refine = normalize_registration(fixed_refine_raw)
moving_refine = normalize_registration(moving_refine_raw)
fixed_mask = tissue_mask(fixed_refine)
moving_mask = tissue_mask(moving_refine)
canvas_center = np.array([(fixed_refine.shape[1] - 1) / 2, (fixed_refine.shape[0] - 1) / 2])


def delta_matrix(parameters: np.ndarray) -> np.ndarray:
    angle, log_sx, log_sy, shear, tx, ty = parameters
    theta = np.deg2rad(angle)
    sx, sy = np.exp(log_sx), np.exp(log_sy)
    linear = np.array(
        [
            [np.cos(theta) * sx, (-np.sin(theta) + shear * np.cos(theta)) * sy, 0],
            [np.sin(theta) * sx, (np.cos(theta) + shear * np.sin(theta)) * sy, 0],
            [0, 0, 1],
        ]
    )
    to_origin = np.array([[1, 0, -canvas_center[0]], [0, 1, -canvas_center[1]], [0, 0, 1]])
    back = np.array([[1, 0, canvas_center[0] + tx], [0, 1, canvas_center[1] + ty], [0, 0, 1]])
    return back @ linear @ to_origin


def refinement_objective(parameters: np.ndarray) -> float:
    total = delta_matrix(parameters) @ initial_refine
    model = transform.AffineTransform(matrix=total)
    warped_mask = transform.warp(
        moving_mask.astype(np.uint8), model.inverse, output_shape=fixed_mask.shape, order=0, preserve_range=True
    ) > 0
    overlap = np.logical_and(fixed_mask, warped_mask)
    if overlap.sum() < 500:
        return 10.0
    warped = transform.warp(moving_refine, model.inverse, output_shape=fixed_refine.shape, order=1, preserve_range=True)
    f, m = fixed_refine[overlap], warped[overlap]
    ncc = float(np.corrcoef(f, m)[0, 1]) if f.std() > 0 and m.std() > 0 else 0.0
    return -(0.78 * dice(fixed_mask, warped_mask) + 0.22 * ncc)


optimized = minimize(
    refinement_objective,
    np.zeros(6),
    method="Powell",
    bounds=[(-4, 4), (-0.06, 0.06), (-0.06, 0.06), (-0.05, 0.05), (-40, 40), (-40, 40)],
    options={"maxiter": 150, "xtol": 0.02, "ftol": 0.0001},
)
print("Registration successful, affine refinement matrix found.", flush=True)
refined_matrix = delta_matrix(optimized.x) @ initial_refine
refined_model = transform.AffineTransform(matrix=refined_matrix)
warped_refine = transform.warp(
    moving_refine, refined_model.inverse, output_shape=fixed_refine.shape, order=1, preserve_range=True
)
warped_refine_mask = transform.warp(
    moving_mask.astype(np.uint8), refined_model.inverse, output_shape=fixed_mask.shape, order=0, preserve_range=True
) > 0
refine_overlap = np.logical_and(fixed_mask, warped_refine_mask)
refine_dice = dice(fixed_mask, warped_refine_mask)
f, m = fixed_refine[refine_overlap], warped_refine[refine_overlap]
refine_ncc = float(np.corrcoef(f, m)[0, 1]) if f.size and f.std() > 0 and m.std() > 0 else 0.0

fixed_dapi, fixed_base_shape = load_channel(FIXED_PATH, OUTPUT_LEVEL, 0)
fixed_sox9, _ = load_channel(FIXED_PATH, OUTPUT_LEVEL, 2)
fixed_krt6b, _ = load_channel(FIXED_PATH, OUTPUT_LEVEL, 3)
moving_dapi, moving_base_shape = load_channel(MOVING_PATH, OUTPUT_LEVEL, 0)
moving_pdgfra, _ = load_channel(MOVING_PATH, OUTPUT_LEVEL, 1)
moving_krt14, _ = load_channel(MOVING_PATH, OUTPUT_LEVEL, 4)

output_matrix = scale_matrix_between_levels(
    refined_matrix,
    moving_refine_raw.shape,
    fixed_refine_raw.shape,
    moving_dapi.shape,
    fixed_dapi.shape,
)
output_model = transform.AffineTransform(matrix=output_matrix)


def warp_u8(array: np.ndarray, order: int = 1) -> np.ndarray:
    warped = transform.warp(
        array,
        output_model.inverse,
        output_shape=fixed_dapi.shape,
        order=order,
        preserve_range=True,
        mode="constant",
        cval=0,
    )
    return np.clip(np.rint(warped), 0, 255).astype(np.uint8)


aligned_dapi = warp_u8(moving_dapi)
aligned_pdgfra = warp_u8(moving_pdgfra)
aligned_krt14 = warp_u8(moving_krt14)
print("Applied affine transform to DAPI, PDGFRa, and KRT14 moving channels.", flush=True)

fixed_mask_output = transform.resize(
    fixed_mask.astype(np.uint8), fixed_dapi.shape, order=0, preserve_range=True, anti_aliasing=False
) > 0
moving_mask_output = transform.warp(
    transform.resize(moving_mask.astype(np.uint8), moving_dapi.shape, order=0, preserve_range=True, anti_aliasing=False),
    output_model.inverse,
    output_shape=fixed_dapi.shape,
    order=0,
    preserve_range=True,
) > 0
common = np.logical_and(fixed_mask_output, moving_mask_output)
yy, xx = np.nonzero(common)
if xx.size == 0:
    def bbox(mask: np.ndarray) -> list[int] | None:
        by, bx = np.nonzero(mask)
        return [int(bx.min()), int(by.min()), int(bx.max()) + 1, int(by.max()) + 1] if bx.size else None
    debug = {
        "optimized_parameters": optimized.x.tolist(),
        "optimized_success": bool(optimized.success),
        "optimized_fun": float(optimized.fun),
        "initial_refine_matrix": initial_refine.tolist(),
        "refined_matrix": refined_matrix.tolist(),
        "output_matrix": output_matrix.tolist(),
        "fixed_mask_bbox": bbox(fixed_mask_output),
        "moving_mask_prewarp_bbox": bbox(
            transform.resize(moving_mask.astype(np.uint8), moving_dapi.shape, order=0, preserve_range=True, anti_aliasing=False) > 0
        ),
        "moving_mask_warped_bbox": bbox(moving_mask_output),
        "fixed_shape": list(fixed_dapi.shape),
        "moving_shape": list(moving_dapi.shape),
    }
    (AUDIT / "empty_overlap_debug.json").write_text(json.dumps(debug, indent=2), encoding="utf-8")
    print(json.dumps(debug, indent=2), flush=True)
    raise RuntimeError("Registered tissue masks do not overlap")
margin = 192
x0 = max(0, int(xx.min()) - margin)
x1 = min(fixed_dapi.shape[1], int(xx.max()) + 1 + margin)
y0 = max(0, int(yy.min()) - margin)
y1 = min(fixed_dapi.shape[0], int(yy.max()) + 1 + margin)
crop = np.s_[y0:y1, x0:x1]

fixed_dapi_n = normalize_display(fixed_dapi[crop])
sox9_n = normalize_display(fixed_sox9[crop])
krt6b_n = normalize_display(fixed_krt6b[crop])
aligned_dapi_n = normalize_display(aligned_dapi[crop])
pdgfra_n = normalize_display(aligned_pdgfra[crop])
krt14_n = normalize_display(aligned_krt14[crop])

lane1 = np.stack(
    [np.clip(krt6b_n.astype(np.uint16) + sox9_n.astype(np.uint16), 0, 255).astype(np.uint8), sox9_n, fixed_dapi_n],
    axis=-1,
)
lane2 = np.stack([krt14_n, pdgfra_n, aligned_dapi_n], axis=-1)
lane3 = np.stack([sox9_n, pdgfra_n, fixed_dapi_n], axis=-1)

fixed_level_scale = fixed_base_shape[1] / fixed_dapi.shape[1]
pixel_size_um = BASE_PIXEL_SIZE_UM * fixed_level_scale
resolution = (10000.0 / pixel_size_um, 10000.0 / pixel_size_um)


def write_tiff(path: Path, array: np.ndarray, rgb: bool) -> None:
    tifffile.imwrite(
        path,
        array,
        bigtiff=True,
        photometric="rgb" if rgb else "minisblack",
        compression="deflate",
        predictor=True,
        tile=(256, 256),
        resolution=resolution,
        resolutionunit="CENTIMETER",
        metadata={"PhysicalSizeX": pixel_size_um, "PhysicalSizeY": pixel_size_um, "PhysicalSizeUnit": "µm"},
    )


write_tiff(ALIGNED / "DAPI_FibK14_aligned_to_KCMig.tif", aligned_dapi[crop], False)
write_tiff(ALIGNED / "PDGFRa_aligned_to_KCMig.tif", aligned_pdgfra[crop], False)
write_tiff(ALIGNED / "KRT14_aligned_to_KCMig.tif", aligned_krt14[crop], False)
write_tiff(LANES / "Lane1_Composite.tif", lane1, True)
write_tiff(LANES / "Lane2_Composite.tif", lane2, True)
write_tiff(LANES / "Lane3_Overlap.tif", lane3, True)
print("Saved Lane1_Composite.tif, Lane2_Composite.tif, and Lane3_Overlap.tif as standard RGB TIFF.", flush=True)


def save_preview(array: np.ndarray, path: Path, max_size: int = 1800) -> None:
    image = Image.fromarray(array)
    image.thumbnail((max_size, max_size), Image.Resampling.LANCZOS)
    image.save(path)


save_preview(lane1, LANES / "Lane1_Composite_preview.png")
save_preview(lane2, LANES / "Lane2_Composite_preview.png")
save_preview(lane3, LANES / "Lane3_Overlap_preview.png")
qc = np.stack(
    [normalize_display(fixed_dapi[crop]), normalize_display(aligned_dapi[crop]), np.zeros_like(fixed_dapi_n)], axis=-1
)
save_preview(qc, AUDIT / "Registration_QC_DAPI_red_fixed_green_aligned.png")

report = {
    "status": "pass" if refine_dice >= 0.75 and refine_ncc >= 0.15 else "review_required",
    "fixed_reference": str(FIXED_PATH),
    "moving": str(MOVING_PATH),
    "true_channel_mapping": {
        "Lane1": "DAPI + SOX9(570) + KRT6B(690); SOX9 rendered yellow",
        "Lane2": "aligned DAPI + PDGFRa(520, green) + KRT14(690, red)",
        "Lane3": "fixed DAPI(blue) + SOX9(570, red) + aligned PDGFRa(520, green)",
    },
    "requested_markers_not_present": ["KRT17", "VIM", "KRT5"],
    "coarse_registration": coarse_report,
    "refinement_level": REFINE_LEVEL,
    "refinement_parameters": optimized.x.tolist(),
    "refined_matrix_moving_to_fixed": refined_matrix.tolist(),
    "output_level": OUTPUT_LEVEL,
    "output_matrix_moving_to_fixed": output_matrix.tolist(),
    "refined_tissue_dice": refine_dice,
    "refined_overlap_ncc": refine_ncc,
    "crop_xyxy_level2": [x0, y0, x1, y1],
    "output_shape_yx": [y1 - y0, x1 - x0],
    "pixel_size_um": pixel_size_um,
    "normalization": "per-channel positive-pixel robust percentile to uint8; no raw files modified",
    "registration_model": "global affine; serial-section differences remain and are not co-localization at the same-cell level",
}
(AUDIT / "registration_and_lane_manifest.json").write_text(
    json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(report, ensure_ascii=False, indent=2))
