#include "streampower.h"

int main(int argc, char** argv)
{
    // default parameter file name
    std::string parameter_file("ThawScape.ini");
    if (argc > 1) {
        // override default parameter file name
        parameter_file = std::string(argv[1]);
    }

	int nx = 10;
	int ny = 10;
	StreamPower sp = StreamPower(nx, ny);
	sp.Init(parameter_file);
	//std::vector<std::vector<float>> topo = sp.CreateRandomField();
	//char* fname = argv[1];
    sp.LoadInputs();
	sp.Start();

}
