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
		void obtain_communities(deque <double>& weak_hist,vector<vector<int> >& quotation,deque<int> tgood);
		void random_obtain_communities(deque <double>& weak_hist,vector<vector<int> >& quotation,deque<int> tgood);
		double process_overlapping(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz);
		double process_overlapping_2(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz);
		double process_overlapping_3(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz);
		double process_overlapping_4(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz_2);

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

void lfm_network::obtain_communities(deque <double>& weak_hist,vector<vector<int> >& quotation,deque<int> tgood){

	set <int> done;
	int done_size=0;
	int numLoopsFind=0;


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
		}
				
	}else{
		nodes_done++;
	}

	

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



double lfm_network::process_overlapping_4(deque <double>& weak_hist,vector<vector<int> >& quotation,vector <int>& big,deque<int>& tbad,vector <double>& fuz_2){
	


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
	
	int k;
	double ob=0;
	double weak_ave=0;
	int t;
	#pragma omp parallel default(none) private(k,t) shared(big,quotation,mem_vector,fuz_2,bigf,weak_hist,tbad) reduction(+:ob,weak_ave)
	{


	
	while (bigf) {
		
		
		#pragma omp barrier
		#pragma omp single
		{
			fuz_2.clear();
			bigf=false;
			fuz_2.assign(quotation.size(),0);
		}
		
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
		{
			
			#pragma omp for schedule(dynamic,32) 
			for ( k=0; k<big.size(); k++) {
			
				for (int j=0; j<quotation[big[k]].size(); j++)
					if (mem_vector[quotation[big[k]][j]]>1){ 
						
							fuz_2[big[k]]++;//sets the number of nodes of each community that belongs to more than one community
					}
			}
			
			
			#pragma omp single
			{
				for (int i=0; i<big.size(); i++) {
				if (!(fuz_2[big[i]]/quotation[big[i]].size()<1)) {
				//neglects the community with overlap > 1 meaning that all the nodes already belong to another community
					
					for (int j=0; j<quotation[big[i]].size(); j++){//size of the community
						mem_vector[ quotation[big[i]] [j] ]--;
						//sets the id of the community to every node of mem_quot
					}
					big.erase(remove(big.begin(),big.end(),big[i]));
					bigf=true; //Sets flag to compute again
					break;
				
					}

				}
			}
			
			}	
		}	
			

		#pragma omp single nowait
		{
		for (t=0; t<dim; t++)
		if (mem_vector[t]==0)//Node doesn't belongs to any community
			tbad.push_back(t);
		}

		#pragma omp for schedule(dynamic,20)
		for (t=0; t<dim; t++) {
			if(mem_vector[t]>1)////Counts the number of nodes that belong to a community
				ob++;
		}
		
		#pragma omp for schedule(dynamic,20)
		for (t=0; t<big.size(); t++)
			weak_ave+=weak_hist[big[t]];

	}

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

	set <int> done;
	int done_size=0;


  	//random_obtain_communities(weak_hist,quotation,tgood);
  	obtain_communities(weak_hist,quotation,tgood);

	int numLoopsFind=0;

	//------------------------------------------------------------------------------------ computes the overlaps; neglects the community with overlap > 1

	
	
	vector <int> big;
	for (int i=0; i<quotation.size(); i++)
		big.push_back(i);

	

	vector <double> fuz_2;

double weak_ave=process_overlapping_4(weak_hist, quotation,big,tbad,fuz_2);


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
			
			
			/*
			if (module_set) {
				for (int y=0; y<v.size(); y++)
					arcout<<vertices[v[y]]->value<<"\t";
				arcout<<-1<<endl;
			}
			
			//*/
			
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




