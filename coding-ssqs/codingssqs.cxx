#include <iostream>                              
#include <fstream>                              
#include <cstdlib>   // exit

using namespace std;

int main( int argc, char* argv[] )
{
	double c_i = 0;
	int      i = 0;
	double a_i, s_i, d_i;
	double d_sum = 0, w_sum = 0, s_sum = 0;

	ifstream data( argc > 1 ? argv[1] : "/dev/null" );
	if( !data ) {
		cerr << "Error opening data file." << endl;
		::exit(1);
	}
	while( data >> a_i >> s_i ) {
		++i;
		if( a_i < c_i ) {
			d_i = c_i - a_i;
			d_sum += d_i;
		} else {
			d_i = 0;
		}
		c_i = a_i + d_i + s_i;

		s_sum += s_i;
		w_sum += d_i + s_i;
	}
	if( i ) {
		cout << "jobs                 " << i << endl;
		cout << "average interarrival " << (a_i/i) << endl;
		cout << "average delay        " << (d_sum/i) << endl;
		cout << "average wait         " << (w_sum/i) << endl;
		cout << "average service      " << (s_sum/i) << endl;
	}
	data.clear();
	data.close();
	return 0;
}

