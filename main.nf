$HOSTNAME = ""
params.outdir = 'results'  



process download_ddbb{
	
	script:
    """
	if [ ! -d ${projectDir}/data ]; then
		mkdir ${projectDir}/data
		cd ${projectDir}/data
		wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/data/ddbb/McPAS-TCR.csv
		wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/data/ddbb/McPAS-TRB_human.csv
		wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/data/ddbb/VDJdb-TCR.tsv
		wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/data/ddbb/VDJdb-TRB_human.tsv
	fi
	"""
}


if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 

Channel.value(params.mate).set{g_42_1_g_61}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_52_0_g_61}
 } else {  
	g_52_0_g_61 = Channel.empty()
 }



process mixcr_analyze {

input:
 set val(name),  file(reads) from g_52_0_g_61
 val mate from g_42_1_g_61

output:
 file "*clonotypes.ALL.txt"  into g_61_outputFileTxt00_g_37
 file "*.report"  into g_61_outputFileTxt11
 file "*.report0[0-4]"  into g_61_outputFileTxt20_g_36
 file "*.clonotypes.*.txt"  into g_61_outputFileTxt30_g_47


script:

readArray = reads.toString().split(' ')	
R1 = readArray.grep(~/.*R1.*/)[0]
R2 = readArray.grep(~/.*R2.*/)[0]



"""
n=\$(grep "${name}[^0-9]" "/work/dolphinnext/vered/run117/sampleslist.csv" | sed 's/,.*//')


mixcr analyze shotgun --species $params.specie --starting-material rna --only-productive \
--align "-OsaveOriginalReads=true" \
${reads} \$n
csplit -f \$n.report \$n.report '/^==/' '{*}' > mixcr_qc.log

if [ ! -d ${projectDir}/files ]; then
  mkdir ${projectDir}/files
fi


cp *clonotypes.ALL.txt ${projectDir}/files
cp *.report ${projectDir}/files
cp *.report0[0-4] ${projectDir}/files
cp *.clonotypes.*.txt ${projectDir}/files
"""
}


process mixcr_qc {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /01_mixcr_qc.html$/) "1_mixcr_qc/$filename"}
input:
 file file from g_61_outputFileTxt20_g_36.collect()

output:
 file "01_mixcr_qc.html"  into g_36_outputFileHTML00
 file "01*{.out,tab,tsv,csv,png}"  into g_36_outputFileTxt10_g_43, g_36_outputFileTxt11_g_58, g_36_outputFileTxt10_g_37


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/01_mixcr_qc.Rmd



Rscript -e "here<-getwd();rmarkdown::render('01_mixcr_qc.Rmd',
    params=list(
        'workDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=FALSE)"

"""

}


process datafiltering {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /02_datafiltering.html$/) "2_datafiltering/$filename"}
input:
 file files_all from g_36_outputFileTxt10_g_37.collect()
 file files_all from g_61_outputFileTxt00_g_37.collect()

output:
 file "02_datafiltering.html"  into g_37_outputFileHTML00
 file "02*{.out,tab,csv,tsv,png}"  into g_37_outputFileTxt10_g_43, g_37_outputFileTxt11_g_58
 file "clones_*_filtered"  into g_37_outputDir21_g_47, g_37_outputDir21_g_53, g_37_outputDir21_g_54, g_37_outputDir21_g_57
 file "clones_*_filtered.Rda"  into g_37_outputFileTxt30_g_48


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/02_datafiltering.Rmd



Rscript -e "here<-getwd();rmarkdown::render('02_datafiltering.Rmd',
    params=list(
        'workDir'=here,
        'inputDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"

"""

}


process dataset_overview {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /03_dataset_overview.html$/) "3_dataset_overview/$filename"}
input:
 file files from g_36_outputFileTxt10_g_43
 file files from g_37_outputFileTxt10_g_43

output:
 file "03_dataset_overview.html"  into g_43_outputFileHTML00
 file "03*{.out,tab,csv,tsv,png}"  into g_43_outputFileTxt10_g_44, g_43_outputFileTxt11_g_58


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/03_dataset_overview.Rmd





Rscript -e "here<-getwd();rmarkdown::render('03_dataset_overview.Rmd',
    params=list(
        'inputDir1'=here,
        'inputDir2'=here,
        'workDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"

"""

}


process correlations {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /04_correlations.html$/) "4_correlations/$filename"}
input:
 file files from g_43_outputFileTxt10_g_44

output:
 file "04_correlations.html"  into g_44_outputFileHTML00
 file "04*{.out,tab,csv,tsv,png}"  into g_44_outputFileTxt10_g_47, g_44_outputFileTxt11_g_58


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/04_correlations.Rmd

 


Rscript -e "here<-getwd();rmarkdown::render('04_correlations.Rmd',
    params=list(
        'workDir'=here,
        'inputDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"

"""

}


process overlap {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /05_overlap.html$/) "5_overlap/$filename"}
input:
 file file from g_44_outputFileTxt10_g_47.collect()
 file file from g_61_outputFileTxt30_g_47.collect()
 file dir from g_37_outputDir21_g_47

output:
 file "05_overlap.html"  into g_47_outputFileHTML00
 file "05*{.out,tab,csv,tsv,png}"  into g_47_outputFileTxt10_g_48, g_47_outputFileTxt11_g_58


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/05_overlap.Rmd




Rscript -e "here<-getwd();rmarkdown::render('05_overlap.Rmd',
    params=list(
        'workDir'=here,
        'inputDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"


"""

}


process diversity {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /06_diversity.html$/) "6_diversity/$filename"}
input:
 file files from g_47_outputFileTxt10_g_48
 file files from g_37_outputFileTxt30_g_48

output:
 file "06_diversity.html"  into g_48_outputFileHTML00
 file "06*{.out,tab,csv,tsv,png}"  into g_48_outputFileTxt10_g_53, g_48_outputFileTxt11_g_58


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/06_diversity.Rmd




Rscript -e "here<-getwd();rmarkdown::render('06_diversity.Rmd',
    params=list(
        'workDir'=here,
        'inputDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"


"""

}


process kmers {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /07_kmers.html$/) "7_kmers/$filename"}
input:
 file files from g_48_outputFileTxt10_g_53
 file dir from g_37_outputDir21_g_53

output:
 file "07_kmers.html"  into g_53_outputFileHTML00
 file "07*{.out,tab,csv,tsv,png}"  into g_53_outputFileTxt10_g_54, g_53_outputFileTxt11_g_58


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/07_kmers.Rmd




Rscript -e "here<-getwd();rmarkdown::render('07_kmers.Rmd',
    params=list(
        'workDir'=here,
        'inputDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"


"""

}


process network {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /08_network.html$/) "8_network/$filename"}
input:
 file files from g_53_outputFileTxt10_g_54
 file dir from g_37_outputDir21_g_54

output:
 file "08_network.html"  into g_54_outputFileHTML00
 file "08*{.out,tab,csv,tsv,png}"  into g_54_outputFileTxt10_g_57, g_54_outputFileTxt11_g_58


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/08_network.Rmd



Rscript -e "here<-getwd();rmarkdown::render('08_network.Rmd',
    params=list(
        'workDir'=here,
        'inputDir'=here,
        'outputDir'=here,
        'sampleInfo'='${projectDir}/sampleslist.csv',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"


"""

}


process ddbb {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /09_ddbb.html$/) "9_ddbb/$filename"}
input:
 file files from g_54_outputFileTxt10_g_57
 file dir from g_37_outputDir21_g_57

output:
 file "09_ddbb.html"  into g_57_outputFileHTML00
 file "09*{.out,tab,csv,tsv,png}"  into g_57_outputFileTxt11_g_58


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/09_ddbb.Rmd



Rscript -e "here<-getwd();rmarkdown::render('09_ddbb.Rmd',
    params=list(
        'workDir'=here,
        'inputDir'=here,
        'outputDir'=here,
        'mcpas'='${projectDir}/${params.mcpas}',
        'vdjdb'='${projectDir}/${params.vdjdb}',
        'sampleInfo'='${projectDir}/sampleslist.csv',
        'chain'='${params.chain}',
        'specie'='${params.specie}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"


"""

}


process report {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /final-report_${params.chain}.html$/) "10_report/$filename"}
input:
 file files from g_36_outputFileTxt11_g_58
 file files from g_43_outputFileTxt11_g_58
 file files from g_44_outputFileTxt11_g_58
 file files from g_48_outputFileTxt11_g_58
 file files from g_53_outputFileTxt11_g_58
 file files from g_54_outputFileTxt11_g_58
 file files from g_57_outputFileTxt11_g_58
 file files from g_47_outputFileTxt11_g_58
 file files from g_37_outputFileTxt11_g_58

output:
 file "final-report_${params.chain}.html"  into g_58_outputFileHTML00


script:

"""
#!/bin/bash

wget https://raw.githubusercontent.com/ConesaLab/TCR_nextflow/main/scripts/10_report.Rmd

Rscript -e "bookdown::render_book(
                input='10_report.Rmd',
                params=list(
                    'inputDir'=getwd(),
                    'chain'='${params.chain}',
                    'specie'='${params.specie}'))"
mv _main.html final-report_${params.chain}.html


"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
