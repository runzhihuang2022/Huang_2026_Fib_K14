import ij.IJ
import java.nio.charset.StandardCharsets


def root = new File('F:/20250325cell背靠背拒稿/20251214上皮化/20260422rebuttal/20260813找曾老师指导/免疫荧光结果统计/Figure 3G overlap/001A_reprocess_20260820/03_imagej_double_positive_spatial')
def auditDir = new File(root, '05_audit')
auditDir.mkdirs()
def relativePaths = [
    '02_recolored_merges/KC_KRT6B_red_SOX9_green_DAPI_blue_RGB.tif',
    '02_recolored_merges/Fib_KRT14_red_PDGFRa_green_DAPI_blue_RGB.tif',
    '03_double_positive_cells/KC_KRT6B_SOX9_double_positive_cell_mask.tif',
    '03_double_positive_cells/Fib_PDGFRa_KRT14_double_positive_cell_mask.tif',
    '03_double_positive_cells/KC_KRT6B_SOX9_double_positive_yellow_DAPI_blue_RGB.tif',
    '03_double_positive_cells/Fib_PDGFRa_KRT14_double_positive_yellow_DAPI_blue_RGB.tif',
    '04_spatial_correspondence/SerialSection_double_positive_cell_centers_RGB.tif',
    '04_spatial_correspondence/SerialSection_10um_cell_neighborhood_RGB.tif',
]
def rows = ['file\tstatus\twidth\theight\tbit_depth\ttype\tpixel_sample']
boolean allPass = true
relativePaths.each { relative ->
    def file = new File(root, relative)
    def imp = IJ.openImage(file.absolutePath)
    if (imp == null) {
        allPass = false
        rows << "${relative}\tFAIL\t\t\t\t\t"
    } else {
        int sample = imp.getProcessor().getPixel(Math.min(100, imp.width - 1), Math.min(100, imp.height - 1))
        String type = imp.bitDepth == 24 ? 'RGB' : (imp.bitDepth == 8 ? 'GRAY8' : "BIT${imp.bitDepth}")
        boolean pass = imp.width == 8544 && imp.height == 3067 && (imp.bitDepth == 8 || imp.bitDepth == 24)
        if (!pass) allPass = false
        rows << "${relative}\t${pass ? 'PASS' : 'FAIL'}\t${imp.width}\t${imp.height}\t${imp.bitDepth}\t${type}\t${sample}"
        imp.close()
    }
}
new File(auditDir, 'fiji_final_output_open_verification.tsv').setText('\uFEFF' + rows.join('\n') + '\n', StandardCharsets.UTF_8.name())
new File(auditDir, 'FIJI_FINAL_OPEN_VERIFIED.txt').setText(allPass ? 'PASS\nfiles=8\n' : 'FAIL\n', StandardCharsets.UTF_8.name())
println(allPass ? 'FINAL_OPEN_VERIFICATION_PASS files=8' : 'FINAL_OPEN_VERIFICATION_FAIL')
System.exit(allPass ? 0 : 1)
