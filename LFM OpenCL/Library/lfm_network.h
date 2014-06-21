#if !defined(LFM_NETWORK)
#define LFM_NETWORK
// DDI class that contains 2 doubles
typedef ddi<int> dd_int;
// iterator of multiset of ddi (double, double)
typedef multiset<dd_int>::iterator fit_it;
//combination of a key value and a mapped value
typedef map <int, fit_it> nodes_map;

	




class lfm_network : public static_network {
	
	
	
	public:
				
		lfm_network(static_network &, int*);
		lfm_network(static_network &, int*,map<int,int> order);
		~lfm_network();
		
		void alpha_setter(double & a, double & b, double & c) {
			initial_alpha=a;
			final_alpha=b;
			alpha_precision=c;
		}
		
		
		int membership(bool);

		
				
	private:
		
		
		// PARAMETERS
		double alpha;
		double initial_alpha;
		double final_alpha;
		double alpha_precision;
		int counter;
		
		double qic;
		double ic_ktot;
		double ic_kin;
		set <int> ic;
		
		double setic(int &);
		double setic_old(int &);
		bool cond (const double &, const double &);
		bool ncond (const double &, const double &);
		double qf(const double &, const double &);	
		
		double set_ics(bool);
		int obtain_communities(deque <double>& weak_hist,vector<vector<int> >& quotation,deque<int> tgood);
		void random_obtain_communities(deque <double>& weak_hist,vector<vector<int> >& quotation,deque<int> tgood);
		double process_overlapping(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz);
		double process_overlapping_2(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz);
		double process_overlapping_opencl(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <int>& fuz_2,int* pos_size_big,int* communities,int num_nodes);

};



bool lfm_network::cond (const double &uno, const double &due) {
	
		
	return ((ic_kin+2.*uno*due)/pow(ic_ktot+due,alpha)-qic>R2_EPS);

}

bool lfm_network::ncond (const double &uno, const double &due) {
			
	return ((ic_kin-2.*uno*due)/pow(ic_ktot-due,alpha)-qic>R2_EPS);

}

double lfm_network::qf(const double &kip, const double &ktp) {
	
		
	return (ic_kin+2.*kip)/pow(ic_ktot+ktp,alpha);

}


lfm_network::lfm_network(static_network &uno, int* seqeuence) : static_network(uno, seqeuence) {

	
	initial_alpha=0.6;
	final_alpha=1.6;
	alpha_precision=0.01;
	
	counter=0;

	
}

lfm_network::lfm_network(static_network &uno, int* seqeuence,map<int,int> order) : static_network(uno, seqeuence,order) {

	
	initial_alpha=0.6;
	final_alpha=1.6;
	alpha_precision=0.01;
	
	counter=0;

	
}



lfm_network::~lfm_network() {}




inline double lfm_network::setic(int &start) {

	qic=0;
	ic_kin=0;
	ic_ktot=0;
	

	nodes_map nodes_fit;
	//multiset multiple elements can have equivalent values
	//dd_init typedef ddi<int> dd_int
	multiset <dd_int> fit_ic;
	multiset <dd_int> fit_con;
	
	//ic is a set
	ic.clear();
	
	int good_node=start;
	//fit_it is an iterator for a multiset
	fit_it good_it= fit_con.insert(dd_int(0., vertices[good_node]->strength, good_node));
	nodes_fit.insert(make_pair(good_node, good_it));

	
	
	
	
	while (good_node!=-1) {
		

		
		
		ic.insert(good_node);
		nodes_fit[good_node]=fit_ic.insert(dd_int(good_it->first, good_it->second, good_node));
		ic_ktot+=vertices[good_node]->strength;		
		ic_kin+=2.*(good_it->first)*(good_it->second);
		qic=(ic_kin)/pow(ic_ktot,alpha);
		fit_con.erase(good_it);

		
		//-------------------------------- updates neighbors
		wsarray* neighbors = vertices[good_node]->links;
		
				
		
		for (int i=0; i<neighbors->size(); i++) {
		
			int vicino=neighbors->l(i);

			nodes_map::iterator it= nodes_fit.find(vicino);
			

			if (it==nodes_fit.end())
				nodes_fit.insert(make_pair(vicino, fit_con.insert(dd_int(neighbors->w(i)/vertices[vicino]->strength, vertices[vicino]->strength, vicino))));
			
			else {
			
				
				if (ic.find(it->first)==ic.end()) {	
				
					double fit_old=(it->second->first)*(it->second->second);
					fit_con.erase(it->second);
					it->second=fit_con.insert(dd_int((neighbors->w(i)+fit_old)/vertices[vicino]->strength, vertices[vicino]->strength, vicino));
				}
				
				
				else {
				
					double fit_old=(it->second->first)*(it->second->second);
					fit_ic.erase(it->second);
					it->second=fit_ic.insert(dd_int((neighbors->w(i)+fit_old)/vertices[vicino]->strength, vertices[vicino]->strength, vicino));
				}

			}
			
		
		}
		
		//--------------------------------  updates neighbors
		
		
		
		//----------------------------------- looks for "bad nodes"
		
		deque <int> bad_nodes;
		fit_it bad_it=fit_ic.begin();
		
		while (ncond(bad_it->first, bad_it->second)) {

			bad_nodes.push_back(bad_it->third);
			bad_it++;

			if (bad_it==fit_ic.end())
				break;
		}
		
		
		//-------------------------------- updates "bad nodes"
		
		
		for (int i=0; i<bad_nodes.size(); i++) {
			
			nodes_map::iterator it__b= nodes_fit.find(bad_nodes[i]);	// it__b->first is the node, it__b->second is the iterator

			if (ncond(it__b->second->first, it__b->second->second)) {
				
				
								
				
				fit_it bad_node_it = it__b->second;
				
				ic.erase(ic.find(bad_nodes[i]));
				it__b->second = fit_con.insert(dd_int(it__b->second->first, it__b->second->second, bad_nodes[i]));
				fit_ic.erase(bad_node_it);
				ic_ktot-=vertices[bad_nodes[i]]->strength;		
				ic_kin-=2.*(it__b->second->first)*(it__b->second->second);
				qic=(ic_kin)/pow(ic_ktot,alpha);

				
				wsarray* neighbors = vertices[bad_nodes[i]]->links;

		
				for (int i=0; i<neighbors->size(); i++) {
				
					int vicino=neighbors->l(i);
					
					nodes_map::iterator it= nodes_fit.find(vicino);
					
					
					if (ic.find(it->first)==ic.end()) {
					
						double fit_old=(it->second->first)*(it->second->second);
						fit_con.erase(it->second);
						it->second=fit_con.insert(dd_int((fit_old-neighbors->w(i))/vertices[vicino]->strength, vertices[vicino]->strength, vicino));
					}
					
					else {	
					
						double fit_old=(it->second->first)*(it->second->second);
						fit_ic.erase(it->second);
						it->second=fit_ic.insert(dd_int((fit_old-neighbors->w(i))/vertices[vicino]->strength, vertices[vicino]->strength, vicino));
						
						if (ncond(it->second->first, it->second->second) && (find(bad_nodes.begin(), bad_nodes.end(), vicino)==bad_nodes.end())) {	// if it is not there
							bad_nodes.push_back(it->first);
							
							
						}
						
					}
				}
			}
				


		
		}
		//-------------------------------- updates bad nodes
		
		
		
		
		

		// ------------- good nodes
		
		
		good_node=-1;

		if (fit_con.empty())
			break;
		
		good_it=fit_con.end();
		good_it--;

		
		
		//------------------------------ randomizes
		
		int choice=0;
		fit_it good_it_2=good_it;
		
		
		deque <fit_it> choices;
		
		while (good_it_2!=fit_con.begin() && (*(good_it_2))==(*(good_it))) {
			
			choice++;
			choices.push_back(good_it_2);
			good_it_2--;
			
		}
		
		if (choices.size()>0)
			good_it=choices[irand(choices.size()-1)];
		
		
		//----------------------------------------------------------------------------------------------------
		
		
		if (cond(good_it->first, good_it->second))
			good_node=good_it->third;
		
			
		
	}
	
	
	
	return (ic_kin/ic_ktot);
	
}



double lfm_network::setic_old(int &start) {
	
		
	set <int> con;
	deque <pair<int, double> > fitcon;

	
	
	deque <int> nps;
	deque <int> bad;
	deque <int> fitic;
	
		
				
	ic_kin=0;
	ic_ktot=0;
	qic=0;
	int nc=0;

	
	
	ic.clear();
	con.insert(start);
	nps.push_back(0);
	pair <int,double> fcp (start, 0);
	fitcon.push_back(fcp);
	fitic.push_back(0);

	
	
	bool flag=true;
	
	while (flag) {
		
		
		
		deque<double> fit(fitcon.size());
		
		
		
		for (int i=0; i<nps.size(); i++) if (i==0 || fit[nps[i]]>0) {
			
			int a=fitcon[nps[i]].first;
			
			
			fitic[nps[i]]=1;
			nc++;

			
			ic_ktot+=vertices[a]->strength;
			ic_kin+= 2* fitcon[nps[i]].second;
			
			qic=qf(0., 0.);

			for (int j=0; j<fitcon.size(); j++)	{					
				
				int b= fitcon[j].first;
				
				
				int indexv=vertices[b]->links->find(a);
				if (indexv!=UNLIKELY)
					fitcon[j].second+=vertices[b]->links->w(indexv);
				
				
				if (fitic[j]==1)
					fit[j]= qic - qf(-fitcon[j].second, -vertices[fitcon[j].first]->strength);

				else
					fit[j]= qf(fitcon[j].second, vertices[fitcon[j].first]->strength) -qic;
					

			}
			
							
				

			for (int j=0; j<vertices[a]->links->size(); j++) {	
							 
				int al=vertices[a]->links->l(j);
				double aw= vertices[a]->links->w(j);
				 
				 
				if (con.insert(al).second) {		
					pair <int,double> fc (al, aw);
					fitcon.push_back(fc);
					fitic.push_back(0);
					fit.push_back(qf(aw, vertices[al]->strength) - qic);
				}
	
			
			}
		}
		
		if (nc==1)
			fit[0]=1;
		
				
		bad.clear();
		nps.clear();
		double npf=0;
		
		for (int i=0; i<fitcon.size(); i++) {
			if (fitic[i]==0) {

				if (fit[i]>npf) {
					nps.clear();
					nps.push_back(i);
					npf=fit[i];
				}
			   
				else if (fit[i]==npf)
					nps.push_back(i);
			}
			
			else if (fit[i]<=0)
				bad.push_back(i);
		
		}
		
		
		if (nps.size()>0) {
		
			int choice=nps[irand(nps.size()-1)];
			nps.clear();
			nps.push_back(choice);
		
		}

		
		for (int i=0; i<bad.size(); i++) if (fit[bad[i]]<=0) {
			
			int a=fitcon[bad[i]].first;
			
			
			fitic[bad[i]]=0;
			nc--;
			ic_ktot-=vertices[a]->strength;
			ic_kin-= 2* fitcon[bad[i]].second;
			
			qic=qf(0., 0.);
			
			for (int j=0; j<fitcon.size(); j++)	{
				
				int b= fitcon[j].first;		
				
				
				int indexv=vertices[b]->links->find(a);
				if (indexv!=UNLIKELY)
					fitcon[j].second-=vertices[b]->links->w(indexv);
				

				if (fitic[j]==1) {			
					fit[j]= qic - qf(-fitcon[j].second, -vertices[fitcon[j].first]->strength);
					
					if (fit[j]<=0 && find(bad.begin(), bad.end(), j)==bad.end())
						bad.push_back(j);
						
									
				}
				else
					fit[j]= qf(fitcon[j].second, vertices[fitcon[j].first]->strength) -qic;
						

			}
				
			
		}
			
		if (bad.size()>0) {
				
			nps.clear();
			npf=0;
				
			for (int i=0; i<fitcon.size(); i++) if(fitic[i]==0){
					
				if (fit[i]>npf) {
					nps.clear();
					nps.push_back(i);
					npf=fit[i];
				}
					
				else if (fit[i]==npf)
					nps.push_back(i);
					
			}
		}
		
		if (nps.size()>0) {
		
			int choice=nps[irand(nps.size()-1)];
			nps.clear();
			nps.push_back(choice);
		
		}
		
		
		
		//-------------------------------------------------------------------------------------------------------------------------------------------------
		
				
		
		flag=(npf>R2_EPS);
		
		
	}
	
	for (int i=0; i<fitic.size(); i++)
		if (fitic[i]==1)
			ic.insert(fitcon[i].first);
	
	return (ic_kin/ic_ktot);
	

	
}

int lfm_network::obtain_communities(deque <double>& weak_hist,vector<vector<int> >& quotation,deque<int> tgood){

	set <int> done;
	int done_size=0;
	int numLoopsFind=0;
	int num_nodes=0;

	
	int nodes_done=0;
	
	//Only if there was an insertion then it means that the node hasn't been
	//assigned to any cover and it finds the community of it
	for (int i=0; i<dim; i++) if(done.insert(i).second) {
		
		ic.clear();
		numLoopsFind++;
		int done_size=done.size();
		double weak_ic;
		// old is always false
		
		weak_ic =setic(i); 
		//obtains the community of the node and
		//saves it on the set ic
			
		//deque <int> ic_id;
		vector <int> ic_id;
		// Iterator of pointer in a set
		set<int>::iterator it = ic.begin();
		//Copy the community form the set to a dequeu
		while(it!=ic.end())
			ic_id.push_back(*(it++));
		
		//Checks if the difference is greater than epsilon, then there was an error
		if (fabs(weak_ic-kin(ic_id)/ktot(ic_id)) > R2_EPS) {
			cout<<"error "<<endl;
			int err;
			cin>>err;
		}

		//only if the community is the same size as the dimension of the network
		if (ic_id.size()==dim)
			tgood.push_back(i);

		//copy all the nodes that have been found		
		for (int w=0; w<ic_id.size(); w++)
			done.insert(ic_id[w]);

		// records the community only if it covers something new
		if (done.size()>done_size) {		
			quotation.push_back(ic_id);
			weak_hist.push_back(weak_ic);
			num_nodes+=ic_id.size();
		}
				
	}else{
		nodes_done++;
	}


	cout<<"Nodes skipped "<<nodes_done<<endl<<endl;
	return num_nodes;
}

void lfm_network::random_obtain_communities(deque <double>& weak_hist,vector<vector<int> >& quotation,deque<int> tgood){

	set <int> done;
	int done_size=0;
	int numLoopsFind=0;
	//Only if there was an insertion then it means that the node hasn't been
	//assigned to any cover and it finds the community of it
	for (int i=0; i<dim; i++) if(done.insert(random_sequence[i]).second) {
		int node=random_sequence[i];
		ic.clear();
		numLoopsFind++;
		int done_size=done.size();
		double weak_ic;
		// old is always false
		
		
			weak_ic =setic(node); 
		//obtains the community of the node and
		//saves it on the set ic
			
		

		
		//deque <int> ic_id;
		vector <int> ic_id;
		// Iterator of pointer in a set
		set<int>::iterator it = ic.begin();
		//Copy the community form the set to a dequeu
		while(it!=ic.end())
			ic_id.push_back(*(it++));
		
		//Checks if the difference is greater than epsilon, then there was an error
		if (fabs(weak_ic-kin(ic_id)/ktot(ic_id)) > R2_EPS) {
			cout<<"error "<<endl;
			int err;
			cin>>err;
		}

		//only if the community is the same size as the dimension of the network
		if (ic_id.size()==dim)
			tgood.push_back(node);

		//copy all the nodes that have been found		
		for (int w=0; w<ic_id.size(); w++)
			done.insert(ic_id[w]);

		// records the community only if it covers something new
		if (done.size()>done_size) {		
			quotation.push_back(ic_id);
			weak_hist.push_back(weak_ic);
		}
				
		
					
	}
}

double lfm_network::process_overlapping(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz_2){
	
	

	deque<deque<int> > mem_quot;			// membership matrix
	vector<int> mem_vector;

	deque <int> first;
	for (int i=0; i<dim; i++)//creates a deque for each vertex
		mem_quot.push_back(first);
	
	vector<int> communities_delete;
	vector<int> nodes_community;

	bool bigf=true;
	
	deque <double> fuz;
	//vector <double> fuz_2;

	
	int times_loop=0;
	while (bigf) {
		
		fuz.clear();
		fuz_2.clear();

		bigf=false;
		int i,k;
		int bigf_int=-1;
		

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		mem_vector.assign(dim,0); 

		fuz_2.assign(quotation.size(),0); 

		//#pragma omp parallel default(none) private (i,k) shared(big,quotation,mem_vector,fuz_2,bigf,bigf_int) 
		{
			//#pragma omp for //default(none) private (i,k) shared(big,quotation,mem_vector)
			for (i=0; i<big.size(); i++) 
				for (int j=0; j<quotation[big[i]].size(); j++){//size of the community
					//#pragma omp atomic
					mem_vector[ quotation[big[i]] [j] ]++;
					//sets the id of the community to every node of mem_quot
				}
			
			//#pragma omp for //default(none) private(i) shared(mem_vector,quotation,big,fuz_2)
			for (k=0; k<big.size(); k++) {
			
				for (int j=0; j<quotation[big[k]].size(); j++)
					if (mem_vector[quotation[big[k]][j]]>1){ 
						
							fuz_2[big[k]]++;//sets the number of nodes of each community that belongs to more than one community
					}
			}
			int counter=0;
			//communities_delete.clear();
			//nodes_community.assign(dim,0);
			//#pragma omp single
			for (int i=0; i<big.size(); i++) {
			if (!(fuz_2[big[i]]/quotation[big[i]].size()<1)) {
			//neglects the community with overlap > 1 meaning that all the nodes already belong to another community
				big.erase(remove(big.begin(),big.end(),big[i]));
				//communities_delete.push_back(big[i]);
				counter++;
				//bigf_int=i;
				bigf=true; //Sets flag to compute again
				break;
			
				}

			}
			
			/*cout<<"Communities Found "<<counter<<endl;

			for (int i=0; i<communities_delete.size(); i++) {
				for (int j=0; j<quotation[communities_delete[i]].size(); j++){
					nodes_community[quotation[ communities_delete[i] ][j]]++;
				}
			}
			int overlap_nodes=0;
			for (int i=0; i<nodes_community.size(); i++) {
				if(nodes_community[i]>0)
					overlap_nodes++;
			}

			cout<<"Overlap Nodes "<<overlap_nodes<<endl;

		
			counter=0;
			for (int i=0; i<big.size(); i++) {
			if (!(fuz_2[big[i]]/quotation[big[i]].size()<1)) {
				//cout<<fuz_2[big[i]]<<"  "<<quotation[big[i]].size()<<endl;
			//neglects the community with overlap > 1 meaning that all the nodes already belong to another community
				big.erase(remove(big.begin(),big.end(),big[i]));
				counter++;
				//bigf=true; //Sets flag to compute again
				//break;
				}
			}

		
		cout<<"Communities Found "<<counter<<endl;*/
		}
			
			
		times_loop++;
			
	}
	

	double ob=0;
	double weak_ave=0;
	int i;
	//#pragma omp parallel default(none) private (i) shared(mem_vector,big,weak_hist,tbad) reduction(+:ob,weak_ave)
	{
		//#pragma omp single nowait
		{
		for (i=0; i<dim; i++)
		if (mem_vector[i]==0)//Node doesn't belongs to any community
			tbad.push_back(i);
		}

		//#pragma omp for
		for (i=0; i<dim; i++) {
			if(mem_vector[i]>1)////Counts the number of nodes that belong to a community
				ob++;
		}
		
		//#pragma omp for
		for (i=0; i<big.size(); i++)
			weak_ave+=weak_hist[big[i]];

	}

weak_ave=weak_ave/big.size();

return weak_ave;

}

double lfm_network::process_overlapping_2(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz_2){
	
	

	deque<deque<int> > mem_quot;			// membership matrix
	vector<int> mem_vector;

	deque <int> first;
	for (int i=0; i<dim; i++)//creates a deque for each vertex
		mem_quot.push_back(first);
	
	vector<int> communities_delete;
	vector<int> nodes_community;

	bool bigf=true;
	
	deque <double> fuz;
	//vector <double> fuz_2;
	
	
	int times_loop=0;

	mem_vector.assign(dim,0); 
	for (int i=0; i<big.size(); i++) 
		for (int j=0; j<quotation[big[i]].size(); j++){//size of the community
			mem_vector[ quotation[big[i]] [j] ]++;
			//sets the id of the community to every node of mem_quot
		}

	while (bigf) {
		
		fuz.clear();
		fuz_2.clear();

		bigf=false;
		int i,k;
		int bigf_int=-1;
		

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		fuz_2.assign(quotation.size(),0); 

		//#pragma omp parallel default(none) private (i,k) shared(big,quotation,mem_vector,fuz_2,bigf,bigf_int) 
		{
			//#pragma omp for //default(none) private (i,k) shared(big,quotation,mem_vector)
			
			
			//#pragma omp for //default(none) private(i) shared(mem_vector,quotation,big,fuz_2)
			for (k=0; k<big.size(); k++) {
			
				for (int j=0; j<quotation[big[k]].size(); j++)
					if (mem_vector[quotation[big[k]][j]]>1){ 
						
							fuz_2[big[k]]++;//sets the number of nodes of each community that belongs to more than one community
					}
			}
			
			
			//#pragma omp single
			for (int i=0; i<big.size(); i++) {
			if (!(fuz_2[big[i]]/quotation[big[i]].size()<1)) {
			//neglects the community with overlap > 1 meaning that all the nodes already belong to another community
				
				for (int j=0; j<quotation[big[i]].size(); j++){//size of the community
					mem_vector[ quotation[big[i]] [j] ]--;
					//sets the id of the community to every node of mem_quot
				}
				big.erase(remove(big.begin(),big.end(),big[i]));
				//communities_delete.push_back(big[i]);
				
				//bigf_int=i;
				bigf=true; //Sets flag to compute again
				break;
			
				}

			}
			
		}	
			
		times_loop++;
			
	}
	

	double ob=0;
	double weak_ave=0;
	int i;
	//#pragma omp parallel default(none) private (i) shared(mem_vector,big,weak_hist,tbad) reduction(+:ob,weak_ave)
	{
		//#pragma omp single nowait
		{
		for (i=0; i<dim; i++)
		if (mem_vector[i]==0)//Node doesn't belongs to any community
			tbad.push_back(i);
		}

		//#pragma omp for
		for (i=0; i<dim; i++) {
			if(mem_vector[i]>1)////Counts the number of nodes that belong to a community
				ob++;
		}
		
		//#pragma omp for
		for (i=0; i<big.size(); i++)
			weak_ave+=weak_hist[big[i]];

	}

weak_ave=weak_ave/big.size();

return weak_ave;

}



double lfm_network::process_overlapping_opencl(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <int>& fuz_2,int* pos_size_big,int* communities,int num_nodes){
	
	vector<int> mem_vector;

	std::vector<cl::Device> devices;     // list of devices found
    cl::Device device;                  // Our chosen device

    cl::Buffer d_communities;                        // device memory used for the input  a vector
    cl::Buffer d_pos_size;                        // device memory used for the input  b vector
    cl::Buffer d_fuz;                       // device memory used for the output c vector
    cl::Buffer d_mem_quot;
    cl::Buffer d_del;

    bool print=false;
	bool bigf=true;
	int times_loop=0;
	int* del=(int*) malloc(sizeof(int)*2);

	del[0]=-1;
	del[1]=-1;

	fuz_2.assign(quotation.size(),1); 
    

	mem_vector.assign(dim,0); 
	for (int i=0; i<quotation.size(); i++) 
		for (int j=0; j<quotation[i].size(); j++){//size of the community
			mem_vector[ quotation[i] [j] ]++;
			//sets the id of the community to every node of mem_quot
		}
	if(print){
		for(int k=0;k<fuz_2.size();k++){
         	cout<<" "<<fuz_2[k];
         }
         cout<<endl<<endl;

	    for(int k=0;k<num_nodes;k++){
	         	cout<<" "<<communities[k];
	         }
	         cout<<endl<<endl;

	     for(int k=0;k<(quotation.size()*2);k++){
	         	cout<<" "<<pos_size_big[k];
	         }
	         cout<<endl<<endl;

		for(int k=0;k<mem_vector.size();k++){
         	cout<<" "<<mem_vector[k];
         }
         cout<<endl<<endl;
       }

	try 
    {
        // Get all platforms
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        // Find a GPU device and set up the context
        cl::Context context;
        for (cl::Platform p : platforms)
        {
            try
            {
                p.getDevices(DEVICE, &devices);
            } catch (cl::Error err) {}// can safely ignore errors here - exception thrown if platform has no GPUs, but >1 platforms
                if (devices.size() > 0)
                {
                    device = devices[0];
                    cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)p(), 0};
                    context = cl::Context(device, properties);
                    break;
                }

        }

         if (!context())
        {
            std::cerr << "Error: Cannot secure a GPU device.\n\nExiting...\n";
            return 1;
        }else{
        	std::cout<<endl<<"GPU secured"<<endl;
        }

         // Load in kernel source, creating a program object for the context

        cl::Program program(context, util::loadProgram("./Library/program.cl"), true);

        // Get the command queue
        cl::CommandQueue queue(context, device);

        auto reset_fuz = cl::make_kernel<cl::Buffer,cl::Buffer>(program, "reset_fuz");

        auto detect_overlap_record = cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer,cl::Buffer>(program, "detect_overlap_record");

        auto erase_community = cl::make_kernel<cl::Buffer,cl::Buffer,cl::Buffer,int,int>(program, "erase_community");

        d_communities   = cl::Buffer(
                    context,
                    CL_MEM_READ_ONLY| CL_MEM_COPY_HOST_PTR,
                    sizeof(int) * num_nodes,
                    &communities[0]);

        d_pos_size   = cl::Buffer(
                    context,
                    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                    sizeof(int) * quotation.size()*2,
                    &pos_size_big[0]);

        d_mem_quot= cl::Buffer(
                    context,
                    CL_MEM_READ_WRITE| CL_MEM_COPY_HOST_PTR,
                    sizeof(int) * dim,
                    &mem_vector[0]);

        d_fuz= cl::Buffer(
                    context,
                    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                    sizeof(int) * quotation.size(),
                    &fuz_2[0]);

        d_del= cl::Buffer(
                    context,
                    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                    sizeof(int) * 2,
                    &del[0]);

         cl::Kernel ko_fuz(program, "reset_fuz");
        ::size_t local = 
            ko_fuz.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);

         std::cout<<"Size of Local "<<local<<endl;
         

         while(bigf){

         	bigf=false;

	         reset_fuz(
	            cl::EnqueueArgs(
	                queue,
	                cl::NDRange(quotation.size()) ), 
	            d_fuz,d_del);


	         detect_overlap_record(
	            cl::EnqueueArgs(
	                queue,
	                cl::NDRange(quotation.size()) ), 
	            d_fuz,
	            d_pos_size,
	            d_communities,
	            d_mem_quot,
	            d_del
	            );
	         
	    	
	    	queue.enqueueReadBuffer(
	            d_del, 
	            CL_TRUE,
	            0, sizeof(int) * 2,
	            &del[0]);

	    	//cout<<"Del Before"<<del[0]<<endl;

	    	if(del[0]!=-1){
	    		//cout<<"Erase Community "<<endl;
	    		int id_com=del[0];
	    		int pos_com=pos_size_big[id_com*2];
	    		int size_com=pos_size_big[(id_com*2) +1];

	    		erase_community(
	            cl::EnqueueArgs(
	                queue,
	                cl::NDRange(size_com) ), 
	            d_pos_size,
	            d_communities,
	            d_mem_quot,
	            id_com,
	            pos_com
	            );

	            del[0]=-1;
	            bigf=true;
	    	}

	    	times_loop++;


    	}

    	queue.enqueueReadBuffer(
            d_fuz, 
            CL_TRUE,
            0, sizeof(int) * quotation.size(),
            &fuz_2[0]);

		queue.enqueueReadBuffer(
            d_pos_size, 
            CL_TRUE,
            0, sizeof(int) * quotation.size()*2,
            &pos_size_big[0]);
		
		queue.enqueueReadBuffer(
            d_mem_quot, 
            CL_TRUE,
            0, sizeof(int) * dim,
            &mem_vector[0]);


    	

    	cout<<"Times Loop "<<times_loop<<endl;
  
    	if(print){
         cout<<"Del After"<<del[0]<<endl;
         cout<<"Fuz "<<endl;
         for(int k=0;k<fuz_2.size();k++){
         	cout<<" "<<fuz_2[k];
         }
         cout<<endl<<endl;
         cout<<"Communities"<<endl;
	    for(int k=0;k<num_nodes;k++){
	         	cout<<" "<<communities[k];
	         }
	         cout<<endl<<endl;
	     cout<<"pos_size_big"<<endl;
	     for(int k=0;k<(quotation.size()*2);k++){
	         	cout<<" "<<pos_size_big[k];
	         }
	         cout<<endl<<endl;
	     cout<<"mem_vector"<<endl;
	     for(int k=0;k<mem_vector.size();k++){
         	cout<<" "<<mem_vector[k];
         }
         cout<<endl<<endl;
     	}
	}
    catch (cl::Error err) {
        std::cout << "Exception\n";
        std::cerr 
            << "ERROR: "
            << err.what()
            << "("
            << err.err()
           << ")"
           << std::endl;


    }

    times_loop=0;
    big.clear();
    //cout<<" Big "<<endl;
    for(int s=0;s<quotation.size();s++){
    	if(pos_size_big[s*2]!=-1){
    		big.push_back(s);
    		//cout<<" "<<s;
    	}
    }
	cout<<endl<<endl;
	cout<<"Size of Big "<<big.size()<<endl;
	/*bigf=true;

	while (bigf) {
		
		
		fuz_2.clear();

		bigf=false;
		int i,k;
		int bigf_int=-1;
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		fuz_2.assign(quotation.size(),0); 

		for (k=0; k<big.size(); k++) {
			
				for (int j=0; j<quotation[big[k]].size(); j++)
					if (mem_vector[quotation[big[k]][j]]>1){ 
						
							fuz_2[big[k]]++;//sets the number of nodes of each community that belongs to more than one community
					}
			}
			
			for (int i=0; i<big.size(); i++) {
			if (( fuz_2[big[i]]==quotation[big[i]].size() )) {
			//if (!((double)fuz_2[big[i]]/quotation[big[i]].size()<1)) {
	
				for (int j=0; j<quotation[big[i]].size(); j++){//size of the community
					mem_vector[ quotation[big[i]] [j] ]--;
					//sets the id of the community to every node of mem_quot
				}
				big.erase(remove(big.begin(),big.end(),big[i]));

				bigf=true; //Sets flag to compute again
				break;
			
				}

			}	
			
		times_loop++;
			
	}

	cout<<"Times Loop "<<times_loop<<endl;
	*/

	double ob=0;
	double weak_ave=0;
	int i;
	
		for (i=0; i<dim; i++)
		if (mem_vector[i]==0)//Node doesn't belongs to any community
			tbad.push_back(i);
		

		
		for (i=0; i<dim; i++) {
			if(mem_vector[i]>1)////Counts the number of nodes that belong to a community
				ob++;
		}
		
		
		for (i=0; i<big.size(); i++)
			weak_ave+=weak_hist[big[i]];

	weak_ave=weak_ave/big.size();

	return weak_ave;

}

double lfm_network::set_ics(bool old) {
	
	// double-ended queue insert and delete in start and end of queue
	//deque doesnt has contiguous memory
	deque <double> weak_hist;	
	//deque<deque<int> > quotation;
	vector<vector<int> > quotation;
	deque<int> tbad;
	deque<int> tgood;

	int* pos_size_big;
	int* communities;

	set <int> done;
	int done_size=0;
	 struct timeval timstr;      /* structure to hold elapsed time */
	  struct rusage ru;           /* structure to hold CPU time--system and user */
	  double tic,toc;             /* floating point numbers to calculate elapsed wallclock time */
	  double usrtim;              /* floating point number to record elapsed user CPU time */
	  double systim;              /* floating point number to record elapsed system CPU time */

	gettimeofday(&timstr,NULL);
  	tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

  	//random_obtain_communities(weak_hist,quotation,tgood);
  	int num_nodes=obtain_communities(weak_hist,quotation,tgood);

	int numLoopsFind=0;
	//------------------------------------------------------------------------------------ computes ic
	
	cout<<"Number loops find communities "<<numLoopsFind<<endl;
	gettimeofday(&timstr,NULL);
	toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	getrusage(RUSAGE_SELF, &ru);
	timstr=ru.ru_utime;        
	usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	timstr=ru.ru_stime;        
	systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	
	char buff[1024];
    sprintf(buff, "Loop ic \n Elapsed time:\t\t\t%.6lf (s)\nElapsed user CPU time:\t\t%.6lf (s)\nElapsed system CPU time:\t%.6lf (s)\n", toc-tic,usrtim,systim);
	std::string time=buff;
	cout<<time<<endl;
	
	gettimeofday(&timstr,NULL);
  	tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

	//------------------------------------------------------------------------------------ computes ic

	communities= (int*)malloc(sizeof(int)*num_nodes);
	pos_size_big=(int*)malloc(sizeof(int)*quotation.size()*2);
	int pos_com=0;
	for (int i=0; i<quotation.size(); i++){
		pos_size_big[i*2]=pos_com;
		pos_size_big[i*2+1]=quotation[i].size();
		for(int g=0;g<quotation[i].size();g++){
			communities[pos_com]=quotation[i][g];
			pos_com++;
		}
	}

	
	vector <int> big;
	vector <int> fuz_2;
	cout<<"Size of big_ "<<quotation.size()<<endl;

	double weak_ave=process_overlapping_opencl(weak_hist, quotation,big,tbad,fuz_2,pos_size_big,communities,num_nodes);

	cout<<"Size of big_ "<<big.size()<<endl;

	gettimeofday(&timstr,NULL);
	toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	getrusage(RUSAGE_SELF, &ru);
	timstr=ru.ru_utime;        
	usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	timstr=ru.ru_stime;
	systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
	
	char buff2[1024];
    sprintf(buff2, "Detect Overlap  \n Elapsed time:\t\t\t%.6lf (s)\nElapsed user CPU time:\t\t%.6lf (s)\nElapsed system CPU time:\t%.6lf (s)\n", toc-tic,usrtim,systim);
	std::string time2=buff2;
	cout<<time2<<endl;

	//	output and archive ----------------------------------------------------------------
	
	
	
	map <double, int>::iterator it_mr = find_value_in_occurrences(weak_ave);
	
	if (weak_ave<1.-R2_EPS) if (it_mr==occurrences.end()) {
		
		archive_pos.insert(make_pair(weak_ave, arcout.tellp()));
		occurrences.insert(make_pair(weak_ave, 1));

		//arcout<<"big modules= "<<big.size()<<"\n--------------------------------------------------------------------\n"<<endl;
		
		for (int i=0; i<big.size(); i++) {
		
			const vector<int> &v= quotation[big[i]];
			
			
			for (int y=0; y<v.size(); y++)
				arcout<<vertices[v[y]]->id_num<<"\t";
			
			arcout<<-1<<endl;
			
			arcout<<"weak = "<<weak_hist[big[i]]<<"\tnumber = "<<i+1<<"\toverlap = "<<fuz_2[big[i]]/quotation[big[i]].size()<<endl;

			arcout<<endl;

		}
		
		arcout<<-2<<endl;

		if (tgood.size()>0) {
		
			
			arcout<<"percolation due to "<<endl;
			const deque<int> &v= tgood;
			
			for (int y=0; y<v.size(); y++)
				arcout<<vertices[v[y]]->id_num<<"\t";
					
			arcout<<endl;
						
		}
		
		if (tbad.size()>0) {
			arcout<<"homeless nodes"<<endl;
			const deque<int> &v= tbad;
			deque <int> module_part;

			for (int y=0; y<v.size(); y++)
				arcout<<vertices[v[y]]->id_num<<"\t";
			
			arcout<<endl;
			
			
			if (module_set) {
				for (int y=0; y<v.size(); y++)
					arcout<<vertices[v[y]]->value<<"\t";
			}
	
			arcout<<endl;

		}
	
		arcout<<"---------------------------------------------"<<endl;
		arcout<<"alpha:\t"<<alpha<<endl;
		arcout<<"number of modules:\t"<<big.size()<<endl;
		//arcout<<"overlap:\t"<<ob<<endl;
		arcout<<"homeless nodes:\t"<<tbad.size()<<endl;
		arcout<<"average weak fitness:\t"<<weak_ave<<endl;
		arcout<<"---------------------------------------------//"<<endl<<endl<<endl;

	}
	
	else {
		
		it_mr->second++;
			
	}
	
	return weak_ave;
	
}



int lfm_network::membership(bool old) {
	
	
	cout<<"alpha set equal to..."<<endl;
	for (alpha=initial_alpha; alpha<final_alpha+R2_EPS; alpha+=alpha_precision) {
		set_ics(old);
		counter++;
		cout<<alpha<<" ";
		
	}
	
	cout<<endl;

	return counter;
}


#endif




