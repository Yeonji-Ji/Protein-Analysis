import pandas as pd
from os.path import basename
from math import sqrt
import argparse

backbone = ["C ", "N ", "CA", "O ", "S "]
filenames = []
pATOM = []
bbATOM = []
lATOM = []
box = {}
InBoxATOM = []
OutBoxATOM = []

class Point:
    #STORES THE X Y AND Z COORDINATES
    def __init__(self, x, y, z):
            self.x = x
            self.y = y
            self.z = z
    #FINDS DISTANCE BETWEEN TWO POINTS
    def Distance(firstpoint, secondPoint):
        d = sqrt(((firstpoint.x-secondPoint.x)**2)+
                ((firstpoint.y-secondPoint.y)**2)+
                ((firstpoint.z-secondPoint.z)**2))
        return d
    def length(self):
        return sqrt((float(self.x)**2)+(float(self.y)**2)+(float(self.z)**2))

class Atom:
  #ATOM CLASS HOLDS INFORMATION READ IN FROM FILE

    def __init__(self,name,idnum,char,tp,att1,att2,x,y,z,av,bv,n):
        self.name = name
        self.idnum = idnum
        self.Element = char
        self.res = tp
        self.chain = att1
        self.resnum = att2
        self.point = Point(float(x),float(y),float(z))
        self.AValue = av
        self.BValue = bv
        self.AtomN = n

def parseProtData(line):
    if "HETATM" in line or "ATOM" in line:
        if not('CL') in line and not('NA') in line and not('OPC') in line and not('T3P') in line and not('ACE') in line and not('NME') in line:
            if str(line[13:15]) in backbone:
        #LINE IS AN INDIVIDUAL LINE IN THE FILE
          #if "ATOM" in line[0:6]:
        # READS IN INFORMATION IF IT IS AN ATOM AND STORES IN ATOM ARRAY
                name = line[0:6]
                idnum = line[6:11]
                char = line[11:16]
                tp = line[16:20]
                att1 = line[20:23]
                att2 =line[23:27]
                x = line[27:38]
                y = line[38:46]
                z = line[46:54]
                av = line[55:60]
                bv = line[60:66]
                n = line[67:79]

                pATOM.append(Atom(name, idnum, char, tp, att1, att2, x, y, z, av, bv, n))
        return
    else:
        return

# def parseBackbone():
#     for atom in pATOM:
#         if atom.Element in backbone:
#             bbATOM.append(atom)

def parseLigData(line):
    if "HETATM" in line or "ATOM" in line:
        if not('CL') in line and not('NA') in line and not('OPC') in line and not('T3P') in line:
        #LINE IS AN INDIVIDUAL LINE IN THE FILE
          #if "ATOM" in line[0:6]:
        # READS IN INFORMATION IF IT IS AN ATOM AND STORES IN ATOM ARRAY
            name = line[0:6]
            idnum = line[6:11]
            char = line[11:16]
            tp = line[16:20]
            att1 = line[20:23]
            att2 =line[23:26]
            x = line[26:38]
            y = line[38:46]
            z = line[46:54]
            av = line[55:61]
            bv = line[60:66]
            n = line[67:79]
            lATOM.append(Atom(name,idnum,char,tp,att1,att2,x,y,z,av,bv,n))
        return
    else:
        return


def boxLigand(size=10):

    x = 0.0
    y = 0.0
    z = 0.0
    maxx = 0.0
    maxy = 0.0
    maxz = 0.0
    minx = 1000.0
    miny = 1000.0
    minz = 1000.0
    count = 0

    for atom in lATOM:
        x+=atom.point.x
        y+=atom.point.y
        z+=atom.point.z
        maxx = max(maxx, atom.point.x)
        minx = min(minx, atom.point.x)
        maxy = max(maxy, atom.point.y)
        miny = min(miny, atom.point.y)
        maxz = max(maxz, atom.point.z)
        minz = min(minz, atom.point.z)
        count+=1

    box["maxxbox"] = maxx + size
    box["minxbox"] = minx - size
    box["maxybox"] = maxy + size
    box["minybox"] = minx - size
    box["maxzbox"] = maxz + size
    box["minzbox"] = minz - size
    return

def findLowBfactor():

    RES_df = pd.DataFrame(columns=["Name", "IDnum", "Element", "Residue", "Chain", "ResNum", "X", "Y", "Z", "AValue", "BValue", "AtomN"])
    maxx = box["maxxbox"]
    minx = box["minxbox"]
    maxy = box["maxybox"]
    miny = box["minybox"]
    maxz = box["maxzbox"]
    minz = box["minzbox"]
    for atom in pATOM:
        if atom.point.x > minx and atom.point.x < maxx and atom.point.y > miny and atom.point.y < maxy and atom.point.z > minz and atom.point.z < maxz:
            InBoxATOM.append(atom)

    for atom in pATOM:
        if atom not in InBoxATOM:
            OutBoxATOM.append(atom)

    for i, atom in enumerate(OutBoxATOM):
        RES_df.loc[i, "Name"] = atom.name
        RES_df.loc[i, "IDnum"] = atom.idnum
        RES_df.loc[i, "Element"] = atom.Element
        RES_df.loc[i, "Residue"] = atom.res
        RES_df.loc[i, "Chain"] = atom.chain
        RES_df.loc[i, "ResNum"] = atom.resnum
        RES_df.loc[i, "X"] = atom.point.x
        RES_df.loc[i, "Y"] = atom.point.y
        RES_df.loc[i, "Z"] = atom.point.z
        RES_df.loc[i, "AValue"] = atom.AValue
        RES_df.loc[i, "BValue"] = atom.BValue
        RES_df.loc[i, "AtomN"] = atom.AtomN
    RES_df = RES_df.sort_values(by="BValue", ascending=True)
    return RES_df

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Find low B factor residues.")
    requirements = parser.add_argument_group("required arguments")
    requirements.add_argument("-p", dest="protein_4amb", help="input file, expected to be in 4amb.pdb format", type=str, required=True)
    requirements.add_argument("-l", dest="ligand", help="input file, expected to be in ligand.pdb format", type=str, required=True)
    parser.add_argument("-d", dest="distance", help="distance from the ligand.", type=int, required=False)
    args = parser.parse_args()

    filenames.append(basename(args.protein_4amb[:-4]))
    filenames.append(basename(args.ligand[:-4]))
    print(args.protein_4amb)
    print(args.ligand)
    print(filenames)

    f = open(args.ligand, "r")
    #GET ALL LINES AND SAVES ALL THE LINES INTO LINES
    lines = f.readlines()
    #LINE COUNT
    for line in lines:
        parseLigData(line)
        if args.distance:
            boxLigand(args.distance)
        else:
            boxLigand()
    f.close()

    f2 = open(args.protein_4amb, "r")
    lines2 = f2.readlines()
    for line in lines2:
        parseProtData(line)
        # parseBackbone()
    df = findLowBfactor()
    if args.distance:
        dist = args.distance
        df.to_csv(filenames[0][:-5]+"-Bfactor-"+str(dist)+"A.csv", ",")
    else:
        df.to_csv(filenames[0][:-5]+"-Bfactor-10A.csv", ",")
    f2.close()