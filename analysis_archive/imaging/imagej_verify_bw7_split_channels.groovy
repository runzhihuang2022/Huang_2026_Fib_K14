import ij.IJ
import java.nio.charset.StandardCharsets

def root = new File('F:/20250325cell背靠背拒稿/20251214上皮化/20260422rebuttal/20260813找曾老师指导/免疫荧光结果统计/Figure 3G overlap/BW7_reprocess_20260820/01_ome_tiff/02_split_channels')
def files = root.listFiles().findAll { it.name.toLowerCase().endsWith('.tif') }.sort { it.name }
def rows = ['file\tstatus\twidth\theight\tbit_depth']
boolean allPass = files.size() == 9
files.each { file ->
    def imp = IJ.openImage(file.absolutePath)
    boolean pass = imp != null && imp.width == 34048 && imp.height == 33978 && imp.bitDepth == 8
    if (!pass) allPass = false
    rows << "${file.name}\t${pass ? 'PASS' : 'FAIL'}\t${imp?.width ?: ''}\t${imp?.height ?: ''}\t${imp?.bitDepth ?: ''}"
    if (imp != null) imp.close()
}
new File(root, 'fiji_split_channel_open_verification.tsv').setText('\uFEFF' + rows.join('\n') + '\n', StandardCharsets.UTF_8.name())
new File(root, 'FIJI_SPLIT_CHANNELS_OPEN_VERIFIED.txt').setText(allPass ? 'PASS\nfiles=9\n' : "FAIL\nfiles=${files.size()}\n", StandardCharsets.UTF_8.name())
println(allPass ? 'SPLIT_CHANNEL_OPEN_PASS files=9' : 'SPLIT_CHANNEL_OPEN_FAIL')
System.exit(allPass ? 0 : 1)
