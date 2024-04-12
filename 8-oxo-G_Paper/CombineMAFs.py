from pathlib import Path

def combine_maf_files(maf_file1, mar_file2):
    maf_path = Path(maf_file1)
    output_file = maf_path.with_name('KBr_combined.maf')
    with open(maf_file1) as m1, open(mar_file2) as m2, open(output_file, 'w') as o:
        header = m1.readline()
        header += m1.readline()
        o.write(header)
        for line in m1:
            o.write(line)
        [m2.readline() for i in range(2)]
        for line in m2:
            o.write(line)


combine_maf_files(  '/media/cam/Data9/CortezAnalysis/Cam_calls/Analysis/Asymmetry/WT/WT_KBr.maf',
                    '/media/cam/Data9/CortezAnalysis/Cam_calls/Analysis/Asymmetry/HMCES/HMCES_KBr.maf')
