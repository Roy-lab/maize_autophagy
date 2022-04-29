import sys

def readName(inname):
	names = []
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		names.append(parts[0])
	f.close()
	return names

def readMod(inname):
	mod = {}
	f = open(inname,'r')
	for l in f:
		parts = l.strip().split('\t')
		mod[parts[0]] = parts[1]
	f.close()
	return mod

def writeMod(outname,mod,names):
	f = open(outname,'w')
	for g in names:
		m = mod.get(g,'-1')
		f.write('%s\t%s\n' % (g,m))
	f.close()

if __name__ == '__main__':
	names = readName(sys.argv[1])
	mod   = readMod(sys.argv[2])
	writeMod(sys.argv[3],mod,names)
