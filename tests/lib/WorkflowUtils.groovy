class WorkflowUtils {

    static String makeManifest(String name, List<String> header, List<Map> rows, String sep="\t" ) {
        def f = new File(name)
        f.parentFile?.mkdirs()

        def headerLine = header.join(sep)
        def lines = rows.collect { row ->
            header.collect { h -> row.get(h, "") }.join(sep)
        }

        f.text = ([headerLine] + lines).join("\n")
        return f.absolutePath
    }

    static String subsetColumns(String filePath, String outPath, int start = 2, int end = 149) {
        def file = new File(filePath)
        def subsetFile = new File(outPath)

        subsetFile.withWriter { writer ->
            file.eachLine { line ->
                def cols = line.split('\t')
                if (cols.size() < end) cols += (cols.size()..end-1).collect { "" }
                writer << cols[(start-1)..(end-1)].join('\t') << "\n"
            }
        }
		return subsetFile.absolutePath
    }
}

