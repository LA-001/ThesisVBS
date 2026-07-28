# Imported from: https://raw.githubusercontent.com/CERN-PH-CMG/cmg-cmssw/CMGTools-from-CMSSW_7_2_3/CMGTools/Production/python/hadd.py
import os
import pprint
import pickle
import shutil

# Check if the specified file is a nanoAOD file.
def checkNano(file):
    from ROOT import TFile
    rf = TFile.Open(file)
    keys = rf.GetListOfKeys()
    if (keys.FindObject("Events") and keys.FindObject("Runs") and keys.FindObject("LuminosityBlocks")):
        isNano = True
    else:
        isNano = False
    rf.Close()
    return isNano


def haddPck(file, odir, idirs):
    sum = None
    for dir in idirs:
        fileName = file.replace( idirs[0], dir )
        pckfile = open(fileName, 'rb')
        obj = pickle.load(pckfile)
        if sum is None:
            sum = obj
        else:
            try:
                sum += obj
            except TypeError:
                pass
                
    oFileName = file.replace( idirs[0], odir )
    pckfile = open(oFileName, 'wb')
    pickle.dump(sum, pckfile)
    txtFileName = oFileName.replace('.pck','.txt')
    txtFile = open(txtFileName, 'w')
    txtFile.write( str(sum) )
    txtFile.write( '\n' )
    txtFile.close()
    

def hadd_batch(ofile, idirs, batch_size=200):
    """Esegue hadd in due passate se ci sono troppi file."""
    import tempfile

    all_files = []
    for dir in idirs:
        f = os.path.join(dir, 'testoutput.root')
        if os.path.exists(f) and os.path.getsize(f) > 1000:
            all_files.append(f)

    if len(all_files) == 0:
        print('WARNING: no valid files found for', ofile)
        return 0

    if len(all_files) <= batch_size:
        cmd = 'hadd -ff -j 8 -n 500 ' + ofile + ' ' + ' '.join(all_files)
        print(cmd)
        return os.system(cmd)

    # Troppi file — merge in gruppi
    print(f'Too many files ({len(all_files)}), merging in batches of {batch_size}...')
    tmpdir = tempfile.mkdtemp()
    partial_files = []

    for i in range(0, len(all_files), batch_size):
        batch = all_files[i:i+batch_size]
        partial = os.path.join(tmpdir, f'partial_{i}.root')
        cmd = 'hadd -ff -j 8 -n 500 ' + partial + ' ' + ' '.join(batch)
        print(cmd)
        exc = os.system(cmd)
        if exc != 0:
            print('---> ABORTING batch <--- with code ' + str(exc))
            shutil.rmtree(tmpdir)
            return exc
        partial_files.append(partial)

    # Merge finale
    cmd = 'hadd -ff -j 8 -n 500 ' + ofile + ' ' + ' '.join(partial_files)
    print(cmd)
    exc = os.system(cmd)
    shutil.rmtree(tmpdir)
    return exc


def haddChunks(idir, removeDestDir=False, cleanUp=False, destdir=None):
    if destdir == None: destdir = idir
        
    chunks = {}
    for file in sorted(os.listdir(idir)):
        filepath = '/'.join([idir, file])
        if os.path.isdir(filepath):
            compdir = file
            try:
                prefix, num = compdir.split('_Chunk')
            except ValueError:
                continue
            chunks.setdefault(prefix, list()).append(filepath)

    if len(chunks) == 0:
        print('warning: no chunk found.')
        return

    for i, (comp, cchunks) in enumerate(chunks.items(), start=1):
        ofile = '/'.join([destdir, comp + '.root'])
        print()
        print("======================")
        print("hadding folder", i, "/", len(chunks))
        print("======================")
        print(ofile, len(cchunks), 'chunks')

        if removeDestDir and os.path.exists(ofile):
            os.remove(ofile)

        exc = hadd_batch(ofile, cchunks)
        if exc != 0:
            print('---> ABORTING <--- with code ' + str(exc))
            exit(1)

    if cleanUp:
        chunkDir = 'Chunks'
        if os.path.isdir('Chunks'):
            shutil.rmtree(chunkDir)
        os.mkdir(chunkDir)
        for comp, cchunks in chunks.items():
            for chunk in cchunks:
                shutil.move(chunk, chunkDir)
        
if __name__ == '__main__':
    import sys
    odir = "/eos/user/l/lallasia/JobOutput"
    haddChunks(odir)
