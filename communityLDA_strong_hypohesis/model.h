/*
 * Copyright (C) 2007 by
 * 
 * 	Xuan-Hieu Phan
 *	hieuxuan@ecei.tohoku.ac.jp or pxhieu@gmail.com
 * 	Graduate School of Information Sciences
 * 	Tohoku University
 *
 * GibbsLDA++ is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * GibbsLDA++ is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

/* 
 * References:
 * + The Java code of Gregor Heinrich (gregor@arbylon.net)
 *   http://www.arbylon.net/projects/LdaGibbsSampler.java
 * + "Parameter estimation for text analysis" by Gregor Heinrich
 *   http://www.arbylon.net/publications/text-est.pdf
 */

#ifndef	_MODEL_H
#define	_MODEL_H

#include "constants.h"
#include "dataset.h"

using namespace std;

// LDA model
class model {
public:
    // fixed options, in set_default_values()
    string wordmapfile;		// file that contains word map [string -> integer id]
    string membermapfile;		// ##add## file that contains members map [string(userid or username) -> integer id]
    string trainlogfile;	// training log file
    string tassignword_suffix;	// ##modify## suffix for topic assignment words
    string tassigncommunity_suffix;	// ##add## suffix for topic assignment communities
    string cassignmember_suffix;	// ##add## suffix for communities assignment members
    string theta_suffix;	// suffix for theta file
    string phi_suffix;		// suffix for phi file
    string lambda_suffix;	// ##add## suffix for lambda file
    string omega_suffix;	// ##add## suffix for omega file
    string others_suffix;	// suffix for file containing other parameters
    string twords_suffix;	// suffix for file containing words-per-topics
    string tcommunities_suffix;	// ##add## suffix for file containing communities-per-topics
    string cmembers_suffix;	// ##add## suffix for file containing members-per-communities
    
    //in parse_args()
    string dir;			// model directory
    string cfile;		// ##modify## content data file
    string ffile;		// ##add## relationship data file
    string model_name;		// model name
    int model_status;		// model status:
				// MODEL_STATUS_UNKNOWN: unknown status
				// MODEL_STATUS_EST: estimating from scratch
				// MODEL_STATUS_ESTC: continue to estimate the model from a previous one
				// MODEL_STATUS_INF: do inference
    
    //in dataset
    dataset * ptrndata;	// pointer to training dataset object
    dataset * pnewdata; // pointer to new dataset object

    mapid2word id2word; // word map [int => string]
    mapid2member id2member; // ##add## member map [int => string]
    
    // --- model parameters and variables ---    
    int U; // ##modify## dataset size (i.e., number of users),in init_est()
    int V; // vocabulary size,in init_est()
    int M; // ##add## members size,in init_est()
    //in parse_args()
    int C; // ##add## number of communities
    int K; // number of topics
    double alpha, beta, gamma, delta; // ##modify## LDA hyperparameters
    int niters; // number of Gibbs sampling iterations
    int liter; // the iteration at which the model was saved
    int savestep; // saving period
    int twords; // print out top words per each topic
    int tcommunities; // print out top communities per each topic
    int cmembers; // print out top members per each community
    int withrawstrs;

    double * p_tw; // ##modify## temp variable for sampling topic-word
    double * p_tc; // ##add## temp variable for sampling topic-community
    double * p_cm; // ##add## temp variable for sampling community-member
    
    //--- latent variables ---
    int ** zw; // ##modify## topic assignments for words, size U x doc.size()
    int ** zc; // ##add## topic assignments for communities, size U x member.size()
    int ** tm; // ##add## communities assignments for members, size U x member.size()
    
    //--- about content(document) ---
    int ** nwt; // ##modify## nwt[i][j]: number of instances of word/term i assigned to topic j, size V x K
    int * nwtsum; // ##modify## nwtsum[j]: total number of words assigned to topic j, size K
    int ** ndt; // ##modify## ndt[i][j]: number of words in document i assigned to topic j, size U x K
    int * ndtsum; // ##modify## ndtsum[i]: total number of words in document i, size U
    
    //--- about community ---
    int ** nct; // ##add## nct[i][j]: number of instances of community i assigned to topic j, size C x K
    int * nctsum; // ##add## nctsum[j]: total number of communities assigned to topic j, size K
    int ** nrt; // ##add## nrt[i][j]: number of communities in the follow_community set i assigned to topic j, size U x K
    int * nrtsum; // ##add## nrtsum[i]: total number of follows in user i, size U, = nfcsum
    
    //--- about member ---
    int ** nmc; // ##add## nmc[i][j]: number of instances of member i assigned to community j, size M x C
    int * nmcsum; // ##add## nmcsum[j]: total number of members assigned to community j, size C
    int ** nrc; // ##add## nrc[i][j]: number of instances of follow_member set i assigned to community j, size U x C
    int * nrcsum; // ##add## nrcsum[i]: total number of follows in user i, size U, = nftsum
    
    //--- parameters ---
    double ** theta; // theta: user-topic distributions, size U x K
    double ** phi; // phi: topic-word distributions, size K x V
    double ** lambda; // ##add## lambda: topic-communitiy distributions, size K x C
    double ** omega; // ##add## omega: community-member distributions, size C x M
    
    // for inference only
    int inf_liter;
    int newU;
    int newV;
    int newM; //##add##
    
    int ** newzw; //##modify##
    int ** newzc; //##add##
    int ** newtm; //##add##
    
    int ** newnwt; // ##modify## 
    int * newnwtsum; // ##modify## 
    int ** newndt; // ##modify## 
    int * newndtsum; // ##modify##
    
    int ** newnct; // ##add## 
    int * newnctsum; // ##add## 
    int ** newnrt; // ##add## 
    int * newnrtsum; // ##add##
    
    int ** newnmc; // ##add## 
    int * newnmcsum; // ##add## 
    int ** newnrc; // ##add## 
    int * newnrcsum; // ##add##
    
    double ** newtheta;
    double ** newphi;
    double ** newlambda;  //##add##
    double ** newomega; //##add##
    // --------------------------------------
    
    model() {
	set_default_values();
    }
          
    ~model();
    
    // set default values for variables
    void set_default_values();   

    // parse command line to get options
    int parse_args(int argc, char ** argv);
    
    // initialize the model
    int init(int argc, char ** argv);
    
    // load LDA model to continue estimating or to do inference
    //int load_model(string model_name);
    
    // save LDA model to files
    // model_name.tassign: topic assignments for words in docs
    // model_name.theta: document-topic distributions
    // model_name.phi: topic-word distributions
    // model_name.others: containing other parameters of the model (alpha, beta, M, V, K)
    int save_model(string model_name);
    int save_model_tassignword(string filename); //##modify##
    int save_model_tassigncommunity(string filename); //##add##
    int save_model_cassignmember(string filename); //##add##
    int save_model_theta(string filename);
    int save_model_phi(string filename);
    int save_model_lambda(string filename); //##add##
    int save_model_omega(string filename); //##add##
    int save_model_others(string filename);
    int save_model_twords(string filename);
    int save_model_tcommunities(string filename); //##add##
    int save_model_cmembers(string filename); //##add##
    /*
    // saving inference outputs
    int save_inf_model(string model_name);
    int save_inf_model_tassignword(string filename); //##modify##
    int save_inf_model_tassigncommunity(string filename); //##add##
    int save_inf_model_cassignmember(string filename); //##add##
    int save_inf_model_newtheta(string filename);
    int save_inf_model_newphi(string filename);
    int save_inf_model_newlambda(string filename); //##add##
    int save_inf_model_newomega(string filename); //##add##
    int save_inf_model_others(string filename);
    int save_inf_model_twords(string filename);
    int save_inf_model_tcommunities(string filename); //##add##
    int save_inf_model_cmembers(string filename); //##add##
    */
    // init for estimation
    int init_est();
    //int init_estc();
	
    // estimate LDA model using Gibbs sampling
    void estimate();
    int sampling_zw(int u, int n); //##modify##
    int sampling_zc(int u, int f); //##add##
    int sampling_tm(int u, int f); //##add##
    void compute_theta();
    void compute_phi();
    void compute_lambda(); //##add##
    void compute_omega(); //##add##
    /*
    // init for inference
    int init_inf();
    // inference for new (unseen) data based on the estimated LDA model
    void inference();
    int inf_sampling_zw(int u, int n); //##modify##
    //int inf_sampling_zc(int u, int m); //##add##
    //int inf_sampling_tm(int u, int m); //##add##
    void compute_newtheta();
    void compute_newphi();
    //void compute_newlambda(); //##add##
    //void compute_newomega(); //##add##
     */
};

#endif

