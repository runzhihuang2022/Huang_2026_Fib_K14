import ij.IJ
import ij.ImagePlus
import ij.gui.Plot
import ij.io.FileSaver
import ij.measure.Measurements
import ij.measure.ResultsTable
import ij.plugin.filter.EDM
import ij.plugin.filter.GaussianBlur
import ij.plugin.filter.ParticleAnalyzer
import ij.process.AutoThresholder
import ij.process.Blitter
import ij.process.ByteProcessor
import ij.process.ColorProcessor
import java.awt.Color
import java.nio.charset.StandardCharsets


def root = new File('F:/20250325cell背靠背拒稿/20251214上皮化/20260422rebuttal/20260813找曾老师指导/免疫荧光结果统计/Figure 3G overlap/BW3_reprocess_20260820/03_imagej_double_positive_spatial')
def inputDir = new File(root, '01_registered_analysis_channels')
def resultDir = new File(root, '03_double_positive_cells')
def spatialDir = new File(root, '04_spatial_correspondence')
def auditDir = new File(root, '05_audit')
[resultDir, spatialDir, auditDir].each { it.mkdirs() }
def logFile = new File(auditDir, 'fiji_cell_analysis.log')
def doneFile = new File(auditDir, 'FIJI_CELL_ANALYSIS_DONE.txt')
def errorFile = new File(auditDir, 'FIJI_CELL_ANALYSIS_ERROR.txt')
logFile.setText('', StandardCharsets.UTF_8.name())
if (doneFile.exists()) doneFile.delete()
if (errorFile.exists()) errorFile.delete()

def log = { String message ->
    println(message)
    logFile.append(message + '\n', StandardCharsets.UTF_8.name())
}

def openByte = { String name ->
    def imp = IJ.openImage(new File(inputDir, name).absolutePath)
    if (imp == null) throw new RuntimeException("Cannot open ${name}")
    def bp = imp.getProcessor().convertToByte(false) as ByteProcessor
    imp.close()
    return bp
}

def crop = { ByteProcessor source, int x, int y, int width, int height ->
    source.setRoi(x, y, width, height)
    return source.crop() as ByteProcessor
}

def subtractSmoothBackground = { ByteProcessor source, double sigma ->
    def work = source.duplicate() as ByteProcessor
    def background = source.duplicate() as ByteProcessor
    new GaussianBlur().blurGaussian(background, sigma, sigma, 0.01)
    work.copyBits(background, 0, 0, Blitter.SUBTRACT)
    return work
}

def diskPercentile = { byte[] pixels, int width, int height, double cx, double cy, int radius, double fraction ->
    int[] hist = new int[256]
    int xMin = Math.max(0, (int)Math.floor(cx - radius))
    int xMax = Math.min(width - 1, (int)Math.ceil(cx + radius))
    int yMin = Math.max(0, (int)Math.floor(cy - radius))
    int yMax = Math.min(height - 1, (int)Math.ceil(cy + radius))
    double r2 = radius * radius
    int count = 0
    for (int y = yMin; y <= yMax; y++) {
        int offset = y * width
        for (int x = xMin; x <= xMax; x++) {
            double dx = x - cx
            double dy = y - cy
            if (dx * dx + dy * dy <= r2) {
                hist[pixels[offset + x] & 0xff]++
                count++
            }
        }
    }
    int target = Math.max(1, (int)Math.ceil(count * fraction))
    int cumulative = 0
    for (int value = 0; value < 256; value++) {
        cumulative += hist[value]
        if (cumulative >= target) return value
    }
    return 255
}

def displayScale = { ByteProcessor source ->
    int[] hist = source.getHistogram()
    long total = 0
    for (int i = 1; i < 256; i++) total += hist[i]
    if (total == 0) return source.duplicate() as ByteProcessor
    def percentile = { double fraction ->
        long target = Math.max(1L, Math.ceil(total * fraction) as long)
        long cumulative = 0
        for (int i = 1; i < 256; i++) {
            cumulative += hist[i]
            if (cumulative >= target) return i
        }
        return 255
    }
    int lo = percentile(0.01) as int
    int hi = percentile(0.998) as int
    if (hi <= lo) hi = Math.min(255, lo + 1)
    def work = source.duplicate() as ByteProcessor
    work.setMinAndMax(lo, hi)
    return work.convertToByte(true) as ByteProcessor
}

def saveGray = { ByteProcessor bp, File file ->
    if (!new FileSaver(new ImagePlus(file.name, bp)).saveAsTiff(file.absolutePath)) throw new RuntimeException("Cannot save ${file}")
}

def saveRgb = { ByteProcessor red, ByteProcessor green, ByteProcessor blue, File file ->
    def cp = new ColorProcessor(red.width, red.height)
    cp.setChannel(1, red)
    cp.setChannel(2, green)
    cp.setChannel(3, blue)
    def imp = new ImagePlus(file.name, cp)
    if (!new FileSaver(imp).saveAsTiff(file.absolutePath)) throw new RuntimeException("Cannot save ${file}")
    def preview = new ImagePlus('preview', cp.duplicate())
    double scale = Math.min(1.0, 1800.0 / Math.max(preview.width, preview.height))
    if (scale < 1.0) IJ.run(preview, 'Size...', "width=${Math.round(preview.width * scale)} height=${Math.round(preview.height * scale)} interpolation=Bilinear")
    new FileSaver(preview).saveAsPng(new File(file.parentFile, file.name.replaceFirst(/\.tif$/, '_preview.png')).absolutePath)
    preview.close()
}

def processPanel = { String panel, ByteProcessor dapi, ByteProcessor markerA, ByteProcessor markerB, String markerAName, String markerBName ->
    double pixelSizeUm = 1.0113780498504639
    int width = dapi.width
    int height = dapi.height
    int coreSize = 1024
    int margin = 32
    int markerRadius = Math.max(2, Math.round(5.0 / pixelSizeUm) as int)
    double backgroundSigma = Math.max(3.0, 15.0 / pixelSizeUm)
    double minAreaPixels = 20.0 / (pixelSizeUm * pixelSizeUm)
    double maxAreaPixels = 300.0 / (pixelSizeUm * pixelSizeUm)
    def cells = []
    int tilesDone = 0
    int tilesTotal = (int)Math.ceil(width / (double)coreSize) * (int)Math.ceil(height / (double)coreSize)
    log("${panel} START size=${width}x${height} tiles=${tilesTotal}")

    for (int coreY = 0; coreY < height; coreY += coreSize) {
        for (int coreX = 0; coreX < width; coreX += coreSize) {
            tilesDone++
            int coreWidth = Math.min(coreSize, width - coreX)
            int coreHeight = Math.min(coreSize, height - coreY)
            int tileX = Math.max(0, coreX - margin)
            int tileY = Math.max(0, coreY - margin)
            int tileEndX = Math.min(width, coreX + coreWidth + margin)
            int tileEndY = Math.min(height, coreY + coreHeight + margin)
            int tileWidth = tileEndX - tileX
            int tileHeight = tileEndY - tileY
            int localCoreX0 = coreX - tileX
            int localCoreY0 = coreY - tileY
            int localCoreX1 = localCoreX0 + coreWidth
            int localCoreY1 = localCoreY0 + coreHeight

            def dapiTile = crop(dapi, tileX, tileY, tileWidth, tileHeight)
            def dapiWork = subtractSmoothBackground(dapiTile, backgroundSigma)
            new GaussianBlur().blurGaussian(dapiWork, 0.8, 0.8, 0.01)
            int dapiThreshold = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, dapiWork.getHistogram())
            dapiWork.setThreshold(Math.min(255, dapiThreshold + 1), 255, 0)
            def nuclearMask = dapiWork.createMask() as ByteProcessor
            new EDM().toWatershed(nuclearMask)
            def rt = new ResultsTable()
            new ParticleAnalyzer(0, Measurements.AREA | Measurements.CENTROID, rt, minAreaPixels, maxAreaPixels).analyze(new ImagePlus('nuclei', nuclearMask))
            if (rt.size() == 0) continue

            def aWork = subtractSmoothBackground(crop(markerA, tileX, tileY, tileWidth, tileHeight), backgroundSigma)
            def bWork = subtractSmoothBackground(crop(markerB, tileX, tileY, tileWidth, tileHeight), backgroundSigma)
            byte[] aPixels = aWork.getPixels() as byte[]
            byte[] bPixels = bWork.getPixels() as byte[]
            for (int row = 0; row < rt.size(); row++) {
                double cx = rt.getValue('X', row)
                double cy = rt.getValue('Y', row)
                if (cx < localCoreX0 || cx >= localCoreX1 || cy < localCoreY0 || cy >= localCoreY1) continue
                int aP90 = diskPercentile(aPixels, tileWidth, tileHeight, cx, cy, markerRadius, 0.90) as int
                int bP90 = diskPercentile(bPixels, tileWidth, tileHeight, cx, cy, markerRadius, 0.90) as int
                cells << [x: cx + tileX, y: cy + tileY, area: rt.getValue('Area', row), a: aP90, b: bP90]
            }
        }
        log("${panel} row=${Math.min(height, coreY + coreSize)}/${height} cells=${cells.size()}")
    }

    int[] histA = new int[256]
    int[] histB = new int[256]
    cells.each { cell -> histA[cell.a as int]++; histB[cell.b as int]++ }
    int thresholdA = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, histA)
    int thresholdB = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, histB)
    cells.each { cell -> cell.doublePositive = (cell.a > thresholdA && cell.b > thresholdB) }
    def doubleCells = cells.findAll { it.doublePositive }
    log("${panel} THRESHOLDS ${markerAName}>${thresholdA} ${markerBName}>${thresholdB}; DAPI=${cells.size()} double_positive=${doubleCells.size()}")

    def centers = new ByteProcessor(width, height)
    def disks = new ByteProcessor(width, height)
    centers.setColor(255)
    disks.setColor(255)
    int radius = Math.max(2, Math.round(5.0 / pixelSizeUm) as int)
    doubleCells.each { cell ->
        int x = Math.max(0, Math.min(width - 1, Math.round(cell.x as double) as int))
        int y = Math.max(0, Math.min(height - 1, Math.round(cell.y as double) as int))
        centers.putPixel(x, y, 255)
        disks.fillOval(x - radius, y - radius, radius * 2 + 1, radius * 2 + 1)
    }
    saveGray(disks, new File(resultDir, "${panel}_double_positive_cell_mask.tif"))
    saveRgb(disks, disks, displayScale(dapi), new File(resultDir, "${panel}_double_positive_yellow_DAPI_blue_RGB.tif"))

    def rows = ["cell_id\tx_px\ty_px\tnuclear_area_px\t${markerAName}_p90\t${markerBName}_p90\t${markerAName}_positive\t${markerBName}_positive\tdouble_positive"]
    cells.eachWithIndex { cell, index ->
        rows << [index + 1, String.format('%.3f', cell.x), String.format('%.3f', cell.y), String.format('%.3f', cell.area), cell.a, cell.b, cell.a > thresholdA, cell.b > thresholdB, cell.doublePositive].join('\t')
    }
    new File(resultDir, "${panel}_per_cell_measurements.tsv").setText('\uFEFF' + rows.join('\n') + '\n', StandardCharsets.UTF_8.name())
    return [panel: panel, cells: cells, doubleCells: doubleCells, centers: centers, disks: disks, thresholdA: thresholdA, thresholdB: thresholdB, markerAName: markerAName, markerBName: markerBName]
}

try {
    log('Fiji/ImageJ DAPI-segmented cell analysis started')
    def kcDapi = openByte('KC_DAPI_fixed.tif')
    def kcSox9 = openByte('KC_SOX9_fixed.tif')
    def kcKrt6b = openByte('KC_KRT6B_fixed.tif')
    def fibDapi = openByte('Fib_DAPI_aligned.tif')
    def fibPdgfra = openByte('Fib_PDGFRa_aligned.tif')
    def fibKrt14 = openByte('Fib_KRT14_aligned.tif')
    def all = [kcDapi, kcSox9, kcKrt6b, fibDapi, fibPdgfra, fibKrt14]
    if (!all.every { it.width == kcDapi.width && it.height == kcDapi.height }) throw new RuntimeException('Channel dimensions do not match')
    log("OPEN_PASS dimensions=${kcDapi.width}x${kcDapi.height}")

    def kc = processPanel('KC_KRT6B_SOX9', kcDapi, kcKrt6b, kcSox9, 'KRT6B', 'SOX9')
    def fib = processPanel('Fib_PDGFRa_KRT14', fibDapi, fibPdgfra, fibKrt14, 'PDGFRa', 'KRT14')

    double pixelSizeUm = 1.0113780498504639
    int proximityPixels = Math.round(10.0 / pixelSizeUm) as int
    def kcNeighborhood = new ByteProcessor(kcDapi.width, kcDapi.height)
    def fibNeighborhood = new ByteProcessor(kcDapi.width, kcDapi.height)
    kcNeighborhood.setColor(255)
    fibNeighborhood.setColor(255)
    kc.doubleCells.each { cell ->
        int x = Math.round(cell.x as double) as int
        int y = Math.round(cell.y as double) as int
        kcNeighborhood.fillOval(x - proximityPixels, y - proximityPixels, proximityPixels * 2 + 1, proximityPixels * 2 + 1)
    }
    fib.doubleCells.each { cell ->
        int x = Math.round(cell.x as double) as int
        int y = Math.round(cell.y as double) as int
        fibNeighborhood.fillOval(x - proximityPixels, y - proximityPixels, proximityPixels * 2 + 1, proximityPixels * 2 + 1)
    }

    def nearestDistancesUm = { List sourceCells, List targetCells ->
        double[] distances = new double[sourceCells.size()]
        for (int i = 0; i < sourceCells.size(); i++) {
            double bestSquared = Double.POSITIVE_INFINITY
            double sx = sourceCells[i].x as double
            double sy = sourceCells[i].y as double
            for (int j = 0; j < targetCells.size(); j++) {
                double dx = sx - (targetCells[j].x as double)
                double dy = sy - (targetCells[j].y as double)
                double squared = dx * dx + dy * dy
                if (squared < bestSquared) bestSquared = squared
            }
            distances[i] = Math.sqrt(bestSquared) * pixelSizeUm
        }
        return distances
    }
    double[] kcNearest = nearestDistancesUm(kc.doubleCells as List, fib.doubleCells as List)
    double[] fibNearest = nearestDistancesUm(fib.doubleCells as List, kc.doubleCells as List)
    double[] radiiUm = [10d, 25d, 50d, 100d, 200d] as double[]
    def percentagesWithin = { double[] distances ->
        double[] values = new double[radiiUm.length]
        for (int i = 0; i < radiiUm.length; i++) {
            int count = 0
            for (double distance : distances) if (distance <= radiiUm[i]) count++
            values[i] = distances.length > 0 ? 100.0 * count / distances.length : Double.NaN
        }
        return values
    }
    double[] kcWithinCurve = percentagesWithin(kcNearest) as double[]
    double[] fibWithinCurve = percentagesWithin(fibNearest) as double[]
    int kcNearFib = 0
    for (double distance : kcNearest) if (distance <= 10.0) kcNearFib++
    int fibNearKc = 0
    for (double distance : fibNearest) if (distance <= 10.0) fibNearKc++
    double kcNearPct = kcNearest.length > 0 ? 100.0 * kcNearFib / kcNearest.length : Double.NaN
    double fibNearPct = fibNearest.length > 0 ? 100.0 * fibNearKc / fibNearest.length : Double.NaN
    def median = { double[] values ->
        if (values.length == 0) return Double.NaN
        double[] sorted = values.clone()
        Arrays.sort(sorted)
        int mid = sorted.length.intdiv(2)
        return sorted.length % 2 == 1 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2.0
    }

    saveRgb(kc.disks as ByteProcessor, fib.disks as ByteProcessor, displayScale(kcDapi), new File(spatialDir, 'SerialSection_double_positive_cell_centers_RGB.tif'))
    saveRgb(kcNeighborhood, fibNeighborhood, displayScale(kcDapi), new File(spatialDir, 'SerialSection_10um_cell_neighborhood_RGB.tif'))

    def plot = new Plot('Registered serial-section double-positive proximity', 'Neighborhood radius (um)', 'Double-positive cells within radius (%)')
    plot.setColor(new Color(210, 35, 45))
    plot.setLineWidth(2)
    plot.addPoints(radiiUm, kcWithinCurve, Plot.LINE)
    plot.addPoints(radiiUm, kcWithinCurve, Plot.CIRCLE)
    plot.setColor(new Color(40, 150, 70))
    plot.addPoints(radiiUm, fibWithinCurve, Plot.LINE)
    plot.addPoints(radiiUm, fibWithinCurve, Plot.CIRCLE)
    plot.setColor(Color.BLACK)
    plot.addLabel(0.05, 0.08, 'Red: KC KRT6B+SOX9+ near Fib cells; Green: Fib PDGFRa+KRT14+ near KC cells')
    plot.setLimits(0, 210, 0, 100)
    new FileSaver(plot.getImagePlus()).saveAsPng(new File(spatialDir, 'double_positive_proximity_radius_curve.png').absolutePath)

    double kcDpPct = kc.cells.size() > 0 ? 100.0 * kc.doubleCells.size() / kc.cells.size() : Double.NaN
    double fibDpPct = fib.cells.size() > 0 ? 100.0 * fib.doubleCells.size() / fib.cells.size() : Double.NaN
    def summary = [
        'panel\tDAPI_nuclei\tdouble_positive_cells\tdouble_positive_percent\tmarker_A_threshold_p90\tmarker_B_threshold_p90',
        "KC_KRT6B_SOX9\t${kc.cells.size()}\t${kc.doubleCells.size()}\t${String.format('%.6f', kcDpPct)}\t${kc.thresholdA}\t${kc.thresholdB}",
        "Fib_PDGFRa_KRT14\t${fib.cells.size()}\t${fib.doubleCells.size()}\t${String.format('%.6f', fibDpPct)}\t${fib.thresholdA}\t${fib.thresholdB}",
    ]
    new File(resultDir, 'double_positive_cell_summary.tsv').setText('\uFEFF' + summary.join('\n') + '\n', StandardCharsets.UTF_8.name())
    def metrics = [
        'metric\tvalue\tinterpretation',
        "KC_double_positive_cells_within_10um_of_Fib_double_positive\t${kcNearFib}/${kc.doubleCells.size()} (${String.format('%.4f', kcNearPct)}%)\tregistered serial-section proximity",
        "Fib_double_positive_cells_within_10um_of_KC_double_positive\t${fibNearKc}/${fib.doubleCells.size()} (${String.format('%.4f', fibNearPct)}%)\tregistered serial-section proximity",
        "KC_to_nearest_Fib_median_distance_um\t${String.format('%.4f', median(kcNearest) as double)}\tdescriptive nearest-neighbor distance",
        "Fib_to_nearest_KC_median_distance_um\t${String.format('%.4f', median(fibNearest) as double)}\tdescriptive nearest-neighbor distance",
    ]
    for (int i = 0; i < radiiUm.length; i++) {
        metrics << "KC_double_positive_within_${radiiUm[i] as int}um_pct\t${String.format('%.4f', kcWithinCurve[i])}\tdescriptive proximity curve"
        metrics << "Fib_double_positive_within_${radiiUm[i] as int}um_pct\t${String.format('%.4f', fibWithinCurve[i])}\tdescriptive proximity curve"
    }
    new File(spatialDir, 'spatial_correspondence_metrics.tsv').setText('\uFEFF' + metrics.join('\n') + '\n', StandardCharsets.UTF_8.name())

    def method = [
        'software=Fiji/ImageJ headless Groovy workflow',
        'pixel_size=1.0113780498504639 um/pixel',
        'DAPI_segmentation=1024-pixel core tiles with 32-pixel overlap; smooth-background subtraction; Gaussian sigma 0.8 pixel; tile-local Otsu; watershed; accepted nuclear area 20-300 um2',
        'marker_measurement=90th percentile after smooth-background subtraction within a 5 um radius disk centered on each accepted DAPI nucleus',
        "KC_locked_thresholds=KRT6B>${kc.thresholdA}; SOX9>${kc.thresholdB}",
        "Fib_locked_thresholds=PDGFRa>${fib.thresholdA}; KRT14>${fib.thresholdB}",
        'double_positive_definition=accepted DAPI nucleus AND both marker p90 values above the panel thresholds',
        'spatial_correspondence=double-positive cell centers in registered adjacent sections within 10 um',
        'evidence_boundary=adjacent serial sections cannot prove same-cell molecular colocalization; report tissue-scale registered spatial proximity/correspondence',
        'threshold_status=exploratory automatic thresholds derived from these two images; confirmatory analysis requires locked thresholds across independent biological replicates',
    ]
    new File(auditDir, 'cell_analysis_method_and_limitations.txt').setText(method.join('\n') + '\n', StandardCharsets.UTF_8.name())
    doneFile.setText("PASS\nKC_DAPI=${kc.cells.size()}\nKC_double_positive=${kc.doubleCells.size()}\nFib_DAPI=${fib.cells.size()}\nFib_double_positive=${fib.doubleCells.size()}\nKC_within_10um_pct=${kcNearPct}\nFib_within_10um_pct=${fibNearPct}\n", StandardCharsets.UTF_8.name())
    log('CELL_ANALYSIS_DONE')
    System.exit(0)
} catch (Throwable error) {
    def sw = new StringWriter()
    error.printStackTrace(new PrintWriter(sw))
    errorFile.setText(sw.toString(), StandardCharsets.UTF_8.name())
    log("CELL_ANALYSIS_FAILED ${error}")
    System.exit(1)
}
