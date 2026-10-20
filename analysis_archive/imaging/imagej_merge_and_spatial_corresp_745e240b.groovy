import ij.IJ
import ij.ImagePlus
import ij.io.FileSaver
import ij.plugin.filter.GaussianBlur
import ij.process.AutoThresholder
import ij.process.Blitter
import ij.process.ByteProcessor
import ij.process.ColorProcessor
import java.nio.charset.StandardCharsets


def root = new File('F:/20250325cell背靠背拒稿/20251214上皮化/20260422rebuttal/20260813找曾老师指导/免疫荧光结果统计/Figure 3G overlap/001A_reprocess_20260820/03_imagej_double_positive_spatial')
def inputDir = new File(root, '01_registered_analysis_channels')
def mergeDir = new File(root, '02_recolored_merges')
def maskDir = new File(root, '03_double_positive_masks')
def spatialDir = new File(root, '04_spatial_correspondence')
def auditDir = new File(root, '05_audit')
[mergeDir, maskDir, spatialDir, auditDir].each { it.mkdirs() }
def logFile = new File(auditDir, 'fiji_execution.log')
def doneFile = new File(auditDir, 'FIJI_DONE.txt')
def errorFile = new File(auditDir, 'FIJI_ERROR.txt')
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

def positivePercentile = { ByteProcessor bp, double fraction ->
    int[] hist = bp.getHistogram()
    long total = 0
    for (int value = 1; value < 256; value++) total += hist[value]
    if (total == 0) return 255
    long target = Math.max(1L, Math.ceil(total * fraction) as long)
    long cumulative = 0
    for (int value = 1; value < 256; value++) {
        cumulative += hist[value]
        if (cumulative >= target) return value
    }
    return 255
}

def displayScale = { ByteProcessor source ->
    int lo = positivePercentile(source, 0.01) as int
    int hi = positivePercentile(source, 0.998) as int
    if (hi <= lo) hi = Math.min(255, lo + 1)
    def work = source.duplicate() as ByteProcessor
    work.setMinAndMax(lo, hi)
    return work.convertToByte(true) as ByteProcessor
}

def backgroundSubtract = { ByteProcessor source, double sigma ->
    def work = source.duplicate() as ByteProcessor
    def background = source.duplicate() as ByteProcessor
    new GaussianBlur().blurGaussian(background, sigma, sigma, 0.01)
    work.copyBits(background, 0, 0, Blitter.SUBTRACT)
    return work
}

def autoMask = { ByteProcessor source, String label ->
    def work = backgroundSubtract(source, 15.0)
    work.smooth()
    int[] hist = work.getHistogram()
    hist[0] = 0
    int threshold = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, hist)
    work.setThreshold(Math.min(255, threshold + 1), 255, 0)
    def mask = work.createMask() as ByteProcessor
    mask.erode()
    mask.dilate()
    log("THRESHOLD ${label} > ${threshold} after Gaussian-background subtraction")
    return [mask: mask, threshold: threshold]
}

def cellMask = { ByteProcessor dapi, String label ->
    def work = backgroundSubtract(dapi, 10.0)
    new GaussianBlur().blurGaussian(work, 0.8, 0.8, 0.01)
    int threshold = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, work.getHistogram())
    work.setThreshold(Math.min(255, threshold + 1), 255, 0)
    def mask = work.createMask() as ByteProcessor
    for (int i = 0; i < 5; i++) mask.dilate()
    log("DAPI_CELL_MASK ${label} > ${threshold}; expanded 5 pixels (~5 um)")
    return [mask: mask, threshold: threshold]
}

def saveGray = { ByteProcessor bp, File file ->
    if (!new FileSaver(new ImagePlus(file.name, bp)).saveAsTiff(file.absolutePath)) {
        throw new RuntimeException("Failed to save ${file}")
    }
}

def saveRgb = { ByteProcessor red, ByteProcessor green, ByteProcessor blue, File file ->
    def cp = new ColorProcessor(red.width, red.height)
    cp.setChannel(1, red)
    cp.setChannel(2, green)
    cp.setChannel(3, blue)
    if (!new FileSaver(new ImagePlus(file.name, cp)).saveAsTiff(file.absolutePath)) {
        throw new RuntimeException("Failed to save ${file}")
    }
    def preview = new ImagePlus('preview', cp.duplicate())
    double scale = Math.min(1.0, 1800.0 / Math.max(preview.width, preview.height))
    if (scale < 1.0) IJ.run(preview, 'Size...', "width=${Math.round(preview.width * scale)} height=${Math.round(preview.height * scale)} interpolation=Bilinear")
    new FileSaver(preview).saveAsPng(new File(file.parentFile, file.name.replaceFirst(/\.tif$/, '_preview.png')).absolutePath)
    preview.close()
}

def countPositive = { ByteProcessor bp ->
    int[] hist = bp.getHistogram()
    long count = 0
    for (int value = 1; value < 256; value++) count += hist[value]
    return count
}

try {
    log('Fiji/ImageJ merge and spatial-correspondence workflow started')
    def kcDapi = openByte('KC_DAPI_fixed.tif')
    def kcSox9 = openByte('KC_SOX9_fixed.tif')
    def kcKrt6b = openByte('KC_KRT6B_fixed.tif')
    def fibDapi = openByte('Fib_DAPI_aligned.tif')
    def fibPdgfra = openByte('Fib_PDGFRa_aligned.tif')
    def fibKrt14 = openByte('Fib_KRT14_aligned.tif')
    def all = [kcDapi, kcSox9, kcKrt6b, fibDapi, fibPdgfra, fibKrt14]
    if (!all.every { it.width == kcDapi.width && it.height == kcDapi.height }) throw new RuntimeException('Channel dimensions do not match')
    log("OPEN_PASS dimensions=${kcDapi.width}x${kcDapi.height}; six analysis channels readable")

    // Fiji recoloring: both panels use red for the keratin, green for the partner marker, blue for DAPI.
    saveRgb(displayScale(kcKrt6b), displayScale(kcSox9), displayScale(kcDapi), new File(mergeDir, 'KC_KRT6B_red_SOX9_green_DAPI_blue_RGB.tif'))
    saveRgb(displayScale(fibKrt14), displayScale(fibPdgfra), displayScale(fibDapi), new File(mergeDir, 'Fib_KRT14_red_PDGFRa_green_DAPI_blue_RGB.tif'))
    log('MERGE_PASS standard RGB TIFFs saved with non-conflicting display colors')

    def kcSox9Result = autoMask(kcSox9, 'KC_SOX9')
    def kcKrt6bResult = autoMask(kcKrt6b, 'KC_KRT6B')
    def kcCellResult = cellMask(kcDapi, 'KC_DAPI')
    def fibPdgfraResult = autoMask(fibPdgfra, 'Fib_PDGFRa')
    def fibKrt14Result = autoMask(fibKrt14, 'Fib_KRT14')
    def fibCellResult = cellMask(fibDapi, 'Fib_DAPI')

    def kcDouble = (kcSox9Result.mask as ByteProcessor).duplicate() as ByteProcessor
    kcDouble.copyBits(kcKrt6bResult.mask as ByteProcessor, 0, 0, Blitter.AND)
    kcDouble.copyBits(kcCellResult.mask as ByteProcessor, 0, 0, Blitter.AND)
    def fibDouble = (fibPdgfraResult.mask as ByteProcessor).duplicate() as ByteProcessor
    fibDouble.copyBits(fibKrt14Result.mask as ByteProcessor, 0, 0, Blitter.AND)
    fibDouble.copyBits(fibCellResult.mask as ByteProcessor, 0, 0, Blitter.AND)
    saveGray(kcDouble, new File(maskDir, 'KC_KRT6B_SOX9_double_positive_mask.tif'))
    saveGray(fibDouble, new File(maskDir, 'Fib_PDGFRa_KRT14_double_positive_mask.tif'))
    saveRgb(kcDouble, kcDouble, displayScale(kcDapi), new File(maskDir, 'KC_double_positive_yellow_DAPI_blue_RGB.tif'))
    saveRgb(fibDouble, fibDouble, displayScale(fibDapi), new File(maskDir, 'Fib_double_positive_yellow_DAPI_blue_RGB.tif'))

    long kcPixels = countPositive(kcDouble) as long
    long fibPixels = countPositive(fibDouble) as long
    def exact = kcDouble.duplicate() as ByteProcessor
    exact.copyBits(fibDouble, 0, 0, Blitter.AND)
    long exactPixels = countPositive(exact) as long
    long unionPixels = kcPixels + fibPixels - exactPixels
    double exactDice = (kcPixels + fibPixels) > 0 ? 2.0 * exactPixels / (kcPixels + fibPixels) : Double.NaN
    double exactJaccard = unionPixels > 0 ? exactPixels / (double)unionPixels : Double.NaN

    def kcDilated10 = kcDouble.duplicate() as ByteProcessor
    def fibDilated10 = fibDouble.duplicate() as ByteProcessor
    for (int i = 0; i < 10; i++) { kcDilated10.dilate(); fibDilated10.dilate() }
    def kcNearFib = kcDouble.duplicate() as ByteProcessor
    kcNearFib.copyBits(fibDilated10, 0, 0, Blitter.AND)
    def fibNearKc = fibDouble.duplicate() as ByteProcessor
    fibNearKc.copyBits(kcDilated10, 0, 0, Blitter.AND)
    long kcNearPixels = countPositive(kcNearFib) as long
    long fibNearPixels = countPositive(fibNearKc) as long
    double kcNearPct = kcPixels > 0 ? 100.0 * kcNearPixels / kcPixels : Double.NaN
    double fibNearPct = fibPixels > 0 ? 100.0 * fibNearPixels / fibPixels : Double.NaN

    saveRgb(kcDouble, fibDouble, displayScale(kcDapi), new File(spatialDir, 'SerialSection_double_positive_spatial_overlay_RGB.tif'))
    saveRgb(kcDilated10, fibDilated10, displayScale(kcDapi), new File(spatialDir, 'SerialSection_10um_neighborhood_overlay_RGB.tif'))

    def rows = [
        'metric\tvalue\tinterpretation',
        "KC_KRT6B_SOX9_double_positive_pixels\t${kcPixels}\texploratory pixel area within DAPI-expanded cell mask",
        "Fib_PDGFRa_KRT14_double_positive_pixels\t${fibPixels}\texploratory pixel area within DAPI-expanded cell mask",
        "exact_overlap_pixels\t${exactPixels}\tregistered serial-section pixel overlap",
        "exact_Dice\t${String.format('%.6f', exactDice)}\tdescriptive only",
        "exact_Jaccard\t${String.format('%.6f', exactJaccard)}\tdescriptive only",
        "KC_pixels_within_10um_of_Fib_double_positive\t${String.format('%.4f', kcNearPct)}%\tspatial proximity after affine registration",
        "Fib_pixels_within_10um_of_KC_double_positive\t${String.format('%.4f', fibNearPct)}%\tspatial proximity after affine registration",
    ]
    new File(spatialDir, 'spatial_correspondence_metrics.tsv').setText('\uFEFF' + rows.join('\n') + '\n', StandardCharsets.UTF_8.name())

    def method = [
        'software=Fiji/ImageJ headless Groovy workflow',
        'color_map_KC=KRT6B red; SOX9 green; DAPI blue',
        'color_map_Fib=KRT14 red; PDGFRa green; DAPI blue',
        'display_scaling=per-channel 1st to 99.8th percentile of positive pixels; display only',
        'marker_mask=Gaussian background subtraction sigma 15 px; smooth; Otsu on nonzero histogram; erode then dilate',
        'cell_support_mask=DAPI Gaussian background subtraction sigma 10 px; Gaussian blur 0.8 px; Otsu; dilate 5 px (~5 um)',
        'double_positive_definition=intersection of both marker masks within the DAPI-expanded cell-support mask',
        'serial_section_proximity=10-pixel (~10 um) neighborhood after the previously audited affine registration',
        'evidence_boundary=adjacent serial sections cannot establish same-cell molecular colocalization; outputs show registered tissue-scale spatial correspondence only',
        'threshold_status=exploratory automatic thresholds; biological review and locked thresholds across independent slides are required for confirmatory statistics',
    ]
    new File(auditDir, 'analysis_method_and_limitations.txt').setText(method.join('\n') + '\n', StandardCharsets.UTF_8.name())
    doneFile.setText("PASS\nwidth=${kcDapi.width}\nheight=${kcDapi.height}\nexact_Dice=${exactDice}\nexact_Jaccard=${exactJaccard}\nKC_within_10um_pct=${kcNearPct}\nFib_within_10um_pct=${fibNearPct}\n", StandardCharsets.UTF_8.name())
    log('ANALYSIS_DONE')
    System.exit(0)
} catch (Throwable error) {
    def sw = new StringWriter()
    error.printStackTrace(new PrintWriter(sw))
    errorFile.setText(sw.toString(), StandardCharsets.UTF_8.name())
    log("ANALYSIS_FAILED ${error}")
    System.exit(1)
}
