
/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *	This program is free software; you can redistribute it and/or modify         *
 *  it under the terms of the GNU General Public License as published by         *
 *  the Free Software Foundation; either version 2 of the License, or            *
 *  (at your option) any later version.                                          *
 *                                                                               *
 *  This program is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
 *  GNU General Public License for more details.                                 *
 *                                                                               *
 *  You should have received a copy of the GNU General Public License            *
 *  along with this program; if not, write to the Free Software                  *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *  Created by Andrea Lancichinetti on 10/07/07 (email: arg.lanci@gmail.com)     *
 *	Modified on 7/07/08                                                          *
 *	Collaborators: Santo Fortunato, Janos Kertesz                                *
 *  Location: ISI foundation, Turin, Italy                                       *
 *	Project: LFM method for community detection                                  *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */




#include "lfm.h"


int write_part(bool, double &, map<int, int> &);
int get_module (double &, ostream &, bool,  map<int, int> &);


int main(int argc, char **argv) {		
	
	
	//-----------------------------------------------------------------------------------------------------
	 struct timeval timstr;      /* structure to hold elapsed time */
	  struct rusage ru;           /* structure to hold CPU time--system and user */
	  double tic,toc;             /* floating point numbers to calculate elapsed wallclock time */
	  double usrtim;              /* floating point number to record elapsed user CPU time */
	  double systim;              /* floating point number to record elapsed system CPU time */
	

	
	srand4();
	
	
	string s;
	char *inputfile = NULL;
	if(argc < 2)
	{
		cerr<<"options: "<<endl;
		cerr<<"giving inputfile using: -f inputfilepath"<<endl;
		exit(0);
	}
	for(int i = 1; i<argc; i++)
	{
		string temp = argv[i];
		if(temp == "-f")
		{
			inputfile = argv[i+1];
			s=inputfile;
			i++;
		}
	}
	/*
	cout<<"Insert input file, please... (karate.dat is an example)"<<endl;
	cin>>s;
	cout<<endl;
	*/

	double initial_alpha=1.0;
	double final_alpha=1.0;
	double alpha_precision=0.1;
	
	int runs=1;

	
	
	
	//-----------------------------------------------------------------------------------------------------
	
	
	bool weighted=false;
	bool value=false;
	
	gettimeofday(&timstr,NULL);
  	tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	
	// reads the input file, and renames the nodes in order to have a sequence from 0 to N-1
	
	if (read_File(weighted, value, s, s)==-1)
		return -1;
	
	//Creates the file format_net.dat
	
	
	{
		//It creates a static_network object with a list of vertices and
		//their connections using format_net.dat
		static_network origin_0("./Library/Files/format_net.dat");
		cout<<"network:: "<<origin_0.size()<<" nodes and "<<origin_0.edges()<<" edges;\t average degree = "<<2*origin_0.edges()/origin_0.size()<<endl;
		//It obtains the largest connected component and writes it into
		//the file connected_component
		origin_0.connected();
	
	}
	
		
	//It creates a static_network object with a list of vertices and
		//their connections using component.dat
	static_network componente("./Library/Files/connected_component.dat");
	cout<<"biggest connected component:: "<<componente.size()<<" nodes and "<<componente.edges()<<" edges ;\t average degree = "<<2*componente.edges()/componente.size()<<endl;
	
	
	componente.traverse(0);



	double total_runs=0;
	
	
	for (int real=0; real<runs; real++) {
	
	
		int sequence[componente.size()];
		// it sets due as a random sequence of integers from 0 to dim-1
		shuffle (sequence, componente.size());
		
		//lfm_network has herecy from static_network
		//Iniitailizes alpha range and precision
		lfm_network ale(componente, sequence,componente.new_order);
		//lfm_network ale(componente, sequence);
		//module_set=f
		ale.set_values(value);
		cout<<"run: "<<initial_alpha<<" "<<final_alpha<<endl;
		//Sets the values and range of alpha
		ale.alpha_setter(initial_alpha, final_alpha, alpha_precision);
		
		bool old=false;
		
		//loop that calls set_ics() using the alpha range
		total_runs+=ale.membership(old);
		cout<<"run: "<<real<<endl;
		
	
		
	}
	
	cout<<"total_runs "<<total_runs<<endl;
	//prints(occurrences);
	
	
	map <int, int> id_value;
	if (value)
		componente.set_id_value(id_value);
	
	write_part(value, total_runs, id_value);
	
	gettimeofday(&timstr,NULL);
	toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	getrusage(RUSAGE_SELF, &ru);
	timstr=ru.ru_utime;        
	usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	timstr=ru.ru_stime;        
	systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	
	char buff[1024];
    sprintf(buff, "Elapsed time:\t\t\t%.6lf (s)\nElapsed user CPU time:\t\t%.6lf (s)\nElapsed system CPU time:\t%.6lf (s)\n", toc-tic,usrtim,systim);
	std::string time=buff;
	cout<<time<<endl;
	
	//*/
	return 0;
	
	
}



int get_module (double & _d_, ostream & pout, bool value, map<int, int> &id_value) {



	map <double, streampos>:: iterator it_pos=archive_pos.find(_d_);
	
	if (it_pos==archive_pos.end()) {
		cerr<<"error in looking for pos"<<endl;
		int err;
		cin>>err;
	}
	
	
	get_part.seekg(it_pos->second);
	
	
	
	int current;
	
	deque<int> deq1;

	while (true) {
		
		
		while (true) {
		
			get_part>>current;

			if (current !=-1)
				deq1.push_back(current);

			else
				break;
		
		
		}
		
		
		sort(deq1.begin(), deq1.end());
		string s;
		getline(get_part, s);
		getline(get_part, s);
		
		
		pout<<s<<endl;
		prints(deq1, pout);
		
		
		if (value) {
			for (int i=0; i<deq1.size(); i++)
				pout<<id_value[deq1[i]]<<"\t";
			pout<<endl;
		}
		
		

		
		pout<<endl;
		
		get_part>>current;
		
		if (current!=-2) {
			
			deq1.clear();
			deq1.push_back(current);
			
		}
		else
			break;

		
	
	
	}


	
	return 0;


}


int write_part(bool value, double & total_runs,  map<int, int> &id_value) {


	double average=0;
	double var=0;
	
	deque < pair <double, double> > ranking;	// normalized occurrence, fitness
	
	
	
	for (map<double, int>:: iterator it=occurrences.begin(); it!=occurrences.end(); it++) {
		
		ranking.push_back(make_pair(	-double(it->second)/total_runs , it->first ));
		average+=double(it->second)/total_runs;
		var+=double(it->second)/total_runs * double(it->second)/total_runs ;
	
	}
	
	average/=occurrences.size();
	var/=occurrences.size();
	
	var-=average*average;
	double dev=sqrt(average);
	
	
	//cout<<"average "<<average << " +/- "<<dev<<endl;
	
	
	sort(ranking.begin(), ranking.end());
	
		
	
	
	//---------------------------------------------------
	
	ofstream fout("output.dat");
	
	fout<<endl<<endl;
	
	fout<<ranking.size()<<" collections of overlapping modules have been detected"<<endl;
	fout<<"their relative occurrences are the following: (decreasing order)"<<endl;
	for (int i=0; i<ranking.size(); i++)
		fout<<i<<"\t"<<-ranking[i].first<<endl;
	
	
	fout<<endl<<endl<<endl;
	
	for (int i=0; i<ranking.size(); i++) {
		
		deque <int> deq;
		string sdeq;
		
		
		fout<<"--------------------------------------------------------------------------------------------------------"<<endl;
		fout<<"collection "<<i<<";\trelative occurrence: "<<-ranking[i].first<<endl;
		get_module(ranking[i].second, fout, value, id_value);
		
				
		string _s_;
		while (_s_!="---------------------------------------------//") {
			getline(get_part, _s_);
			fout<<_s_<<endl;
		}
		
		fout<<endl<<endl<<endl;
	
	
	}
	
	
	
	
	
	return 0;

}

