cwlVersion: v1.0
class: CommandLineTool

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: RNA-Seq_analysis
      entry:
        $include: ./RNA-Seq_analysis.R

baseCommand: ["Rscript",--save,RNA-Seq_analysis]

inputs: []

outputs: 
  output_file:
    type: File
    outputBinding:
      glob: DEGs.xlsx