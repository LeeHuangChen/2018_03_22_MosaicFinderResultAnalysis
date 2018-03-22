import Configurations as conf
import os
from cPickle import load
from src import util


def readResultFile(folder, filename):
    fusedDict = {}
    familyDict = {}

    with open(os.path.join(folder, filename)) as f:
        for i, line in enumerate(f):
            if i > 0:
                # reference:
                # 0         1               2               3                   4
                # gene_id   fusion_event    family_number   start_med_align1    end_med_align1
                # 5                 6               7
                # start_med_align2  end_med_align2  break_point

                arr = line.split("\t")

                # family proteins only have 3 values
                pid = arr[0]
                fEvent = arr[1]
                familyNumMF = arr[2]
                if pid[0] == 'f':
                    # print line
                    splitInfo = pid.split("_")
                    if len(splitInfo) > 3:
                        familyNumber = -2
                    else:
                        familyNumber = int(splitInfo[0][1:])

                    if fEvent in familyDict.keys():
                        subDict = familyDict[fEvent]
                        if type(subDict[1]) != list:
                            subDict[1] = []
                            subDict[2] = []
                            subDict[3] = []
                        subDict[int(familyNumMF)].append(familyNumber)
                    else:
                        familyDict[fEvent] = {0: [], 1: [], 2: [], 3: []}
                        familyDict[fEvent][int(familyNumMF)].append(familyNumber)
                else:
                    # print line, len(arr)
                    splitInfo = pid.split("_")
                    if len(splitInfo) <= 2:
                        f1 = -1
                        f2 = -1
                    else:
                        # id = splitInfo[1]
                        f1 = int(splitInfo[3])
                        f2 = int(splitInfo[5])
                        # g = int(splitInfo[7])
                        # split = int(splitInfo[9])
                    fusedDict[pid] = [fEvent, f1, f2]
    return fusedDict, familyDict

def main():
    filenames = os.listdir(conf.inputFolder)
    util.generateDirectories(conf.resultFolder)
    with open(os.path.join(conf.resultFolder, "Results.txt"), "w") as wf:
        wf.write("TEvo\tNFam\tNFusions\tavgConf\tnumProt\n")
        for filename in filenames:
            # reference:
            # 0 1    2    3    4    5 6        7 8    9   10   11 12
            # M_mjtt_SeqL_1000_NFam_2_NFusions_2_TEvo_1.5_NGen_5_ BorderInformation
            parsed = filename.split("_")
            # model = parsed[1]
            # seqLen = parsed[3]
            NFam = parsed[5]
            NFusions = parsed[7]
            TEvo = parsed[9]
            # NGen = parsed[11]

            fusedDict, familyDict = readResultFile(conf.inputFolder, filename)

            confidenceArr = []

            for pid in fusedDict.keys():
                fEvent, f1, f2 = fusedDict[pid]

                mf1s = familyDict[fEvent][1]
                mf2s = familyDict[fEvent][2]
                success = 0
                totalAssigns = len(mf1s) + len(mf2s)
                print pid
                print mf1s
                print mf2s
                #print familyDict[fEvent][3]
                success1 = 0
                success2 = 0
                for f in mf1s:
                    if f == f1:
                        success1 += 1
                for f in mf2s:
                    if f == f1:
                        success2 += 1
                if success1 > success2:
                    mf1 = mf1s
                    mf2 = mf2s
                else:
                    mf1 = mf2s
                    mf2 = mf1s

                for f in mf1:
                    if f == f1:
                        success1 += 1
                for f in mf2:
                    if f == f2:
                        success2 += 1

                # if totalAssigns == 0:
                #     confi = 0
                # else:
                #     confi = float(success)/totalAssigns

                if success1 > 0 and success2 > 0:
                    confi = 1
                elif success1 == 0 and success2 == 0:
                    confi = 0
                else:
                    confi = .5

                confidenceArr.append(confi)
                print confi, success, totalAssigns

            avgConf = reduce(lambda x, y: x + y, confidenceArr)/float(len(confidenceArr))

            wf.write(str(TEvo) + "\t" + str(NFam) + "\t" + str(NFusions) + "\t")
            wf.write(str(avgConf) + "\t" + str(len(fusedDict.keys())) + "\n")


if __name__ == '__main__':
    main()
