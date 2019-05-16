#include "streampower.h"

int main(int argc, char** argv)
{
	int nx = 10;
	int ny = 10;
	StreamPower sp = StreamPower(nx, ny);
	sp.Init();
	//std::vector<std::vector<float>> topo = sp.CreateRandomField();
	//char* fname = argv[1];
	sp.SetTopo(sp.ReadArcInfoASCIIGrid("topo.asc"));
	sp.SetFA(sp.ReadArcInfoASCIIGrid("FA.asc"));
	//sp.SetFA(sp.ReadArcInfoASCIIGrid("SedThickness.asc"));   // Option to set sediment thickness
	sp.Start();

}