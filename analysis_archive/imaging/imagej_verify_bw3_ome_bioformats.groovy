import loci.formats.ImageReader
import loci.formats.FormatTools
import java.nio.charset.StandardCharsets

def root = new File('F:/20250325cell背靠背拒稿/20251214上皮化/20260422rebuttal/20260813找曾老师指导/免疫荧光结果统计/Figure 3G overlap/BW3_reprocess_20260820/01_ome_tiff')
def audit = new File(root, '00_conversion_audit')
audit.mkdirs()
def targets = [
    [name: 'Fib_K14_BW3 PDGFRα 520+P63 570+KLF4 620+KRT14 690(1).ome.tiff', channels: 5, sx: 38144, sy: 37914],
    [name: 'KC_Mig_浅ⅡBW3.ome.tiff', channels: 4, sx: 41984, sy: 41901],
]
def rows = ['file\tstatus\tseries\tsize_x\tsize_y\tsize_c\tpixel_type\tbase_plane_probe_bytes']
boolean allPass = true
targets.each { spec ->
    def reader = new ImageReader()
    try {
        def path = new File(root, spec.name as String).absolutePath
        reader.setId(path)
        int seriesCount = reader.getSeriesCount()
        reader.setSeries(0)
        int sx = reader.getSizeX(), sy = reader.getSizeY(), sc = reader.getSizeC()
        int w = Math.min(64, sx), h = Math.min(64, sy)
        byte[] probe = reader.openBytes(0, 0, 0, w, h)
        String pt = FormatTools.getPixelTypeString(reader.getPixelType())
        boolean pass = seriesCount >= 1 && sx == (spec.sx as int) && sy == (spec.sy as int) && sc == (spec.channels as int) && probe.length > 0
        if (!pass) allPass = false
        rows << "${spec.name}\t${pass ? 'PASS' : 'FAIL'}\t${seriesCount}\t${sx}\t${sy}\t${sc}\t${pt}\t${probe.length}"
    } catch (Throwable e) {
        allPass = false
        rows << "${spec.name}\tFAIL\t\t\t\t\t\t${e.class.simpleName}: ${e.message}"
    } finally {
        try { reader.close() } catch (Throwable ignored) {}
    }
}
new File(audit, 'fiji_bioformats_open_verification.tsv').setText('\uFEFF' + rows.join('\n') + '\n', StandardCharsets.UTF_8.name())
new File(audit, 'FIJI_BIOFORMATS_OPEN_VERIFIED.txt').setText(allPass ? 'PASS\nfiles=2\n' : 'FAIL\n', StandardCharsets.UTF_8.name())
println(allPass ? 'BIOFORMATS_OPEN_VERIFICATION_PASS files=2' : 'BIOFORMATS_OPEN_VERIFICATION_FAIL')
System.exit(allPass ? 0 : 1)
