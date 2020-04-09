cwlVersion: v1.1
class: CommandLineTool
inputs:
  metadata: File[]
  metadataSchema: File
  subjects: string[]
outputs:
  merged: stdout
stdout: mergedmetadata.ttl
requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: {$include: merge-metadata.py}
        entryname: merge-metadata.py
baseCommand: [python3, merge-metadata.py]
