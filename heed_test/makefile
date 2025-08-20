default:
	cd ./build && make && ./sim

all:
	cd ./build && make && ./sim && ./shower

shower:
	cd ./build && make && ./shower

quick:
	cd ./build && make && ./quick

hole:
	cd ./build && make && ./gem_hole


root:
	# root -e 'TFile *f = TFile::Open("./build/data.root"); clusters->Draw("z:y:x", "", "");'
	root -e 'new TBrowser'
