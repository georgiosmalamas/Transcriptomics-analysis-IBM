cwlVersion: v1.0
class: CommandLineTool

requirements:
  InitialWorkDirRequirement:
    listing:
    - entryname: test_microarrays_analysis
      entry:
        $include: ./test_microarrays_analysis.R

baseCommand: ["Rscript",--save,test_microarrays_analysis]

inputs: []

outputs: 
  output_file:
    type: File
    outputBinding:
      glob: DEGs.xlsx