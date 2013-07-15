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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"
#include "model.h"

using namespace std;

model::~model() {
    if (p_tw) {
	delete p_tw;
    }
    
    if (p_tc) {
        delete p_tc;
    }
    
    if (p_cm) {
        delete p_cm;
    }

    if (ptrndata) {
	delete ptrndata;
    }
    
    if (pnewdata) {
	delete pnewdata;
    }
    
    //--- latent variables ---
    if (zw) {
	for (int u = 0; u < U; u++) {
	    if (zw[u]) {
		delete zw[u];
	    }
	}
    }
    
    if (zc) {
        for (int u = 0; u < U; u++) {
            if (zc[u]) {
                delete zc[u];
            }
        }
    }
    
    if (tm) {
        for (int u = 0; u < U; u++) {
            if (tm[u]) {
                delete tm[u];
            }
        }
    }
    
    //--- about content(document) ---
    if (nwt) {
	for (int w = 0; w < V; w++) {
	    if (nwt[w]) {
		delete nwt[w];
	    }
	}
    }
    
    if (nwtsum) {
        delete nwtsum;
    }
    
    if (ndt) {
	for (int u = 0; u < U; u++) {
	    if (ndt[u]) {
		delete ndt[u];
	    }
	}
    }
    
    if (ndtsum) {
        delete ndtsum;
    }
    
    //--- about community ---
    if (nct) {
        for (int c = 0; c < C; c++) {
            if (nct[c]) {
                delete nct[c];
            }
        }
    }
    
    if (nctsum) {
        delete nctsum;
    }
    
    if (nrt) {
        for (int u = 0; u < U; u++) {
            if (nrt[u]) {
                delete nrt[u];
            }
        }
    }
    
    if (nrtsum) {
        delete nrtsum;
    }
    
    //--- about member ---
    if (nmc) {
        for (int m = 0; m < M; m++) {
            if (nmc[m]) {
                delete nmc[m];
            }
        }
    }
    
    if (nmcsum) {
        delete nmcsum;
    }
    
    if (nrc) {
        for (int u = 0; u < U; u++) {
            if (nrc[u]) {
                delete nrc[u];
            }
        }
    }
    
    if (nrcsum) {
        delete nrcsum;
    }
    
    //--- parameters ---
    if (theta) {
	for (int u = 0; u < U; u++) {
	    if (theta[u]) {
		delete theta[u];
	    }
	}
    }
    
    if (phi) {
	for (int k = 0; k < K; k++) {
	    if (phi[k]) {
		delete phi[k];
	    }
	}
    }
    
    if (lambda) {
        for (int k = 0; k < K; k++) {
            if (lambda[k]) {
                delete lambda[k];
            }
        }
    }
    
    if (omega) {
        for (int c = 0; c < C; c++) {
            if (omega[c]) {
                delete omega[c];
            }
        }
    }

    // only for inference
    if (newzw) {
	for (int u = 0; u < newU; u++) {
	    if (newzw[u]) {
            delete newzw[u];
	    }
	}
    }
    
    if (newzc) {
        for (int u = 0; u < newU; u++) {
            if (newzc[u]) {
                delete newzc[u];
            }
        }
    }
    
    if (newtm) {
        for (int u = 0; u < newU; u++) {
            if (newtm[u]) {
                delete newtm[u];
            }
        }
    }
    
    //--- about content(document) ---
    if (newnwt) {
        for (int w = 0; w < newV; w++) {
            if (newnwt[w]) {
                delete newnwt[w];
            }
        }
    }
    
    if (newnwtsum) {
        delete newnwtsum;
    }
    
    if (newndt) {
        for (int u = 0; u < newU; u++) {
            if (newndt[u]) {
                delete newndt[u];
            }
        }
    }
    
    if (newndtsum) {
        delete newndtsum;
    }
    
    //--- about community ---
    if (newnct) {
        for (int c = 0; c < C; c++) {
            if (newnct[c]) {
                delete newnct[c];
            }
        }
    }
    
    if (newnctsum) {
        delete newnctsum;
    }
    
    if (newnrt) {
        for (int u = 0; u < newU; u++) {
            if (newnrt[u]) {
                delete newnrt[u];
            }
        }
    }
    
    if (newnrtsum) {
        delete newnrtsum;
    }
    
    //--- about member ---
    if (newnmc) {
        for (int m = 0; m < newM; m++) {
            if (newnmc[m]) {
                delete newnmc[m];
            }
        }
    }
    
    if (newnmcsum) {
        delete newnmcsum;
    }
    
    if (newnrc) {
        for (int u = 0; u < newU; u++) {
            if (newnrc[u]) {
                delete newnrc[u];
            }
        }
    }
    
    if (newnrcsum) {
        delete newnrcsum;
    }
    
    //--- parameters ---
    if (newtheta) {
        for (int u = 0; u < newU; u++) {
            if (newtheta[u]) {
                delete newtheta[u];
            }
        }
    }
    
    if (newphi) {
        for (int k = 0; k < K; k++) {
            if (newphi[k]) {
                delete newphi[k];
            }
        }
    }
    
    if (newlambda) {
        for (int k = 0; k < K; k++) {
            if (newlambda[k]) {
                delete newlambda[k];
            }
        }
    }
    
    if (newomega) {
        for (int c = 0; c < C; c++) {
            if (newomega[c]) {
                delete newomega[c];
            }
        }
    }
}

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
    membermapfile = "membermap.txt";
    trainlogfile = "trainlog.txt";
    tassignword_suffix = ".tassignword";
    tassigncommunity_suffix = ".tassigncommunity";
    cassignmember_suffix = ".cassignmember";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    lambda_suffix = ".lambda";
    omega_suffix = ".omega";
    others_suffix = ".others";
    twords_suffix = ".twords";
    tcommunities_suffix = ".tcommunities";
    cmembers_suffix = ".cmembers";
    
    dir = "./";
    cfile = "content.dat";
    ffile = "follow.dat";
    model_name = "model-final";    
    model_status = MODEL_STATUS_UNKNOWN;
    
    ptrndata = NULL;
    pnewdata = NULL;
    
    U = 0;
    V = 0;
    M = 0;
    C = 100;
    K = 100;
    alpha = 50.0 / K;
    beta = gamma = delta = 0.1;
    niters = 2000;
    liter = 0;
    savestep = 200;    
    twords = 0;
    tcommunities = 0;
    cmembers = 0;
    withrawstrs = 0;
    
    p_tw = NULL;
    p_tc = NULL;
    p_cm = NULL;
    zw = NULL;
    zc = NULL;
    tm = NULL;
    nwt = NULL;
    nwtsum = NULL;
    ndt = NULL;
    ndtsum = NULL;
    nct = NULL;
    nctsum = NULL;
    nrt = NULL;
    nrtsum = NULL;
    nmc = NULL;
    nmcsum = NULL;
    nrc = NULL;
    nrcsum =NULL;
    theta = NULL;
    phi = NULL;
    lambda = NULL;
    omega = NULL;
    
    newU = 0;
    newV = 0;
    newM = 0;
    newzw = NULL;
    newzc = NULL;
    newtm = NULL;
    newnwt = NULL;
    newnwtsum = NULL;
    newndt = NULL;
    newndtsum = NULL;
    newnct = NULL;
    newnctsum = NULL;
    newnrt = NULL;
    newnrtsum = NULL;
    newnmc = NULL;
    newnmcsum = NULL;
    newnrc = NULL;
    newnrcsum =NULL;
    newtheta = NULL;
    newphi = NULL;
    newlambda = NULL;
    newomega = NULL;
}

int model::parse_args(int argc, char ** argv) {
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {
    // call parse_args
    if (parse_args(argc, argv)) {
	return 1;
    }
    
    if (model_status == MODEL_STATUS_EST) {
	// estimating the model from scratch
	if (init_est()) {
	    return 1;
	}
	
    } /*else if (model_status == MODEL_STATUS_ESTC) {
	// estimating the model from a previously estimated one
	if (init_estc()) {
	    return 1;
	}
	
    } else if (model_status == MODEL_STATUS_INF) {
	// do inference
	if (init_inf()) {
	    return 1;
	}
    }*/
   
    else {
        printf("Now the tool is only privide fuction -est");
    }
    return 0;
}

int model::init_est() {
    int u, n, f, m, w, k, c;  //U, Nu, Fu, M, V, K, C
    
    p_tw = new double[K];   //temp variable for sampling topic-word 
    p_tc = new double[K];   //temp variable for sampling topic-community
    p_cm = new double[C];   //temp variable for sampling community-member
    
    // + read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + cfile, dir + ffile, dir + wordmapfile, dir + membermapfile)) {
        printf("Fail to read training data!\n");
        return 1;
    }
    
    // + allocate memory and assign values for variables
    U = ptrndata->U;
    V = ptrndata->V;
    M = ptrndata->M;
    
    // K and C: from command line or default value
    // alpha, beta, gamma, delta: from command line or default values
    // niters, savestep: from command line or default values
    
    nwt = new int*[V];
    for (w = 0; w < V; w++) {
        nwt[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nwt[w][k] = 0;
        }
    }
	
    nwtsum = new int[K];
    for (k = 0; k < K; k++) {
        nwtsum[k] = 0;
    }
    
    ndt = new int*[U];
    for (u = 0; u < U; u++) {
        ndt[u] = new int[K];
        for (k = 0; k < K; k++) {
    	    ndt[u][k] = 0;
        }
    }
    
    ndtsum = new int[U];
    for (u = 0; u < U; u++) {
        ndtsum[u] = 0;
    }
    
    nct = new int*[C];
    for (c = 0; c < C; c++) {
        nct[c] = new int[K];
        for (k = 0; k < K; k++) {
    	    nct[c][k] = 0;
        }
    }
	
    nctsum = new int[K];
    for (k = 0; k < K; k++) {
        nctsum[k] = 0;
    }
    
    nrt = new int*[U];
    for (u = 0; u < U; u++) {
        nrt[u] = new int[K];
        for (k = 0; k < K; k++) {
    	    nrt[u][k] = 0;
        }
    }
    
    nrtsum = new int[U];
    for (u = 0; u < U; u++) {
        nrtsum[u] = 0;
    }
    
    nmc = new int*[M];
    for (m = 0; m < M; m++) {
        nmc[m] = new int[C];
        for (c = 0; c < C; c++) {
    	    nmc[m][c] = 0;
        }
    }
	
    nmcsum = new int[C];
    for (c = 0; c < C; c++) {
        nmcsum[c] = 0;
    }
    
    nrc = new int*[U];
    for (u = 0; u < U; u++) {
        nrc[u] = new int[C];
        for (c = 0; c < C; c++) {
    	    nrc[u][c] = 0;
        }
    }
    
    nrcsum = new int[U];
    for (u = 0; u < U; u++) {
        nrcsum[u] = 0;
    }
    
    srandom(time(0)); // initialize for random number generation
    zw = new int*[U];
    zc = new int*[U];
    tm = new int*[U];
    for (u = 0; u < ptrndata->U; u++) {
        int N = ptrndata->users[u]->doc->length;
        int F = ptrndata->users[u]->fol->length;
        zw[u] = new int[N];
        zc[u] = new int[F];
        tm[u] = new int[F];
        
        // initialize for zw
        for (n = 0; n < N; n++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    zw[u][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    nwt[ptrndata->users[u]->doc->words[n]][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    ndt[u][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwtsum[topic] += 1;
        }
        // total number of words in document i
        ndtsum[u] = N;
        
        // initialize for tm, must be before zc
        for (f = 0; f < F; f++) {
    	    int community = (int)(((double)random() / RAND_MAX) * C);
            
            ptrndata->users[u]->fol->communities[f] = community;
    	    tm[u][f] = community;
    	    
    	    // number of instances of member i assigned to community j
    	    nmc[ptrndata->users[u]->fol->members[f]][community] += 1;
    	    // number of members in user i assigned to community j
    	    nrc[u][community] += 1;
    	    // total number of members assigned to community j
    	    nmcsum[community] += 1;
        }
        // total number of follows in users i
        nrcsum[u] = F;
        
        // initialize for zc
        for (f = 0; f < F; f++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    zc[u][f] = topic;
    	    
    	    // number of instances of community i assigned to topic j
    	    nct[ptrndata->users[u]->fol->communities[f]][topic] += 1;
    	    // number of communities in user i assigned to topic j
    	    nrt[u][topic] += 1;
    	    // total number of communities assigned to topic j
    	    nctsum[topic] += 1;
        }
        // total number of follows in users i
        nrtsum[u] = F;
        
    }
    
    theta = new double*[U];
    for (u = 0; u < U; u++) {
        theta[u] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }
    
    lambda = new double*[K];
    for (k = 0; k < K; k++) {
        lambda[k] = new double[C];
    }
    
    omega = new double*[C];
    for (c = 0; c < C; c++) {
        omega[c] = new double[M];
    }
    
    return 0;
}
/*
int model::init_estc() {
    // estimating the model from a previously estimated one
    int m, n, w, k;
    
    p = new double[K];
    
    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
        printf("Fail to load word-topic assignmetn file of the model!\n");
        return 1;
    }
    
    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
        nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
        ndsum[m] = 0;
    }
    
    for (m = 0; m < ptrndata->M; m++) {
        int N = ptrndata->docs[m]->length;
        
        // assign values for nw, nd, nwsum, and ndsum
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        }
        // total number of words in document i
        ndsum[m] = N;
    }
	
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    
    
    return 0;        
}
*/
/*
int model::init_inf() {
    // estimating the model from a previously estimated one
    int m, n, w, k;
    
    p = new double[K];
    
    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
        printf("Fail to load word-topic assignmetn file of the model!\n");
        return 1;
    }
    
    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
        nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
        ndsum[m] = 0;
    }
    
    for (m = 0; m < ptrndata->M; m++) {
        int N = ptrndata->docs[m]->length;
        
        // assign values for nw, nd, nwsum, and ndsum
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        }
        // total number of words in document i
        ndsum[m] = N;
    }
    
    // read new data for inference
    pnewdata = new dataset;
    if (withrawstrs) {
        if (pnewdata->read_newdata_withrawstrs(dir + cfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
        }
    } else {
        if (pnewdata->read_newdata(dir + cfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
        }
    }
    
    newM = pnewdata->M;
    newV = pnewdata->V;
    
    newnw = new int*[newV];
    for (w = 0; w < newV; w++) {
        newnw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnw[w][k] = 0;
        }
    }
	
    newnd = new int*[newM];
    for (m = 0; m < newM; m++) {
        newnd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnd[m][k] = 0;
        }
    }
	
    newnwsum = new int[K];
    for (k = 0; k < K; k++) {
        newnwsum[k] = 0;
    }
    
    newndsum = new int[newM];
    for (m = 0; m < newM; m++) {
        newndsum[m] = 0;
    }
    
    srandom(time(0)); // initialize for random number generation
    newz = new int*[newM];
    for (m = 0; m < pnewdata->M; m++) {
        int N = pnewdata->docs[m]->length;
        newz[m] = new int[N];
        
        // assign values for nw, nd, nwsum, and ndsum
        for (n = 0; n < N; n++) {
    	    int w = pnewdata->docs[m]->words[n];
    	    int _w = pnewdata->_docs[m]->words[n];
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    newz[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    newnw[_w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    newnd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    newnwsum[topic] += 1;
        }
        // total number of words in document i
        newndsum[m] = N;
    }
    
    newtheta = new double*[newM];
    for (m = 0; m < newM; m++) {
        newtheta[m] = new double[K];
    }
	
    newphi = new double*[K];
    for (k = 0; k < K; k++) {
        newphi[k] = new double[newV];
    }    
    
    return 0;        
}
*/

void model::estimate() {
    if (twords > 0) {
        // print out top words per topic
        dataset::read_wordmap(dir + wordmapfile, &id2word);
    }
    if (cmembers > 0) {
        // print out top communities per topic
        dataset::read_membermap(dir + membermapfile, &id2member);
    }
    
    printf("Sampling %d iterations!\n", niters);
    
    int last_iter = liter;
    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
        printf("Iteration %d ...\n", liter);
        
        // for all zw_i
        for (int u = 0; u < U; u++) {
            for (int n = 0; n < ptrndata->users[u]->doc->length; n++) {
                // (zw_i = zw[m][n])
                // sample from p(zw_i|zw_-i, w)
                int topic = sampling_zw(u, n);
                zw[u][n] = topic;
            }
            
            for (int f = 0; f < ptrndata->users[u]->fol->length; f++) {
                // sample from p(zc_i|zc_-i, c)
                int topic = sampling_zc(u, f);
                zc[u][f] = topic;
            }
            
            for (int f = 0; f < ptrndata->users[u]->fol->length; f++) {
                // sample from p(zc_i|zc_-i, c)
                int community = sampling_tm(u, f);
                tm[u][f] = community;
            }
        }
        
        if (savestep > 0) {
            if (liter % savestep == 0) {
                // saving the model
                printf("Saving the model at iteration %d ...\n", liter);
                compute_theta();
                compute_phi();
                compute_lambda();
                compute_omega();
                save_model(utils::generate_model_name(liter));
            }
        }
    }
    
    printf("Gibbs sampling completed!\n");
    printf("Saving the final model!\n");
    compute_theta();
    compute_phi();
    compute_lambda();
    compute_omega();
    liter--;
    save_model(utils::generate_model_name(-1));
}

int model::sampling_zw(int u, int n) {
    // remove zw_i from the count variables
    int topic = zw[u][n];
    int w = ptrndata->users[u]->doc->words[n];
    nwt[w][topic] -= 1;//number of word assigned to topic
    ndt[u][topic] -= 1;//number of words in doucument(user's content) assigned to topic
    nwtsum[topic] -= 1;//total words in topic
    ndtsum[u] -= 1;//total words in document(users' content)
    
    double Vbeta = V * beta;
    double Kalpha = K * alpha;
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
        p_tw[k] = (nwt[w][k] + beta) / (nwtsum[k] + Vbeta) *
        (ndt[u][k] + nrt[u][k] + alpha) / (ndtsum[u] + nrtsum[u] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
        p_tw[k] += p_tw[k - 1];
    }
    // scaled sample because of unnormalized p_tw[]
    double u_value = ((double)random() / RAND_MAX) * p_tw[K - 1];
    
    for (topic = 0; topic < K; topic++) {
        if (p_tw[topic] > u_value) {
            break;
        }
    }
    
    // add newly estimated zw_i to count variables
    nwt[w][topic] += 1;
    ndt[u][topic] += 1;
    nwtsum[topic] += 1;
    ndtsum[u] += 1;
    
    return topic;
}

int model::sampling_zc(int u, int f) {
    // remove zc_i from the count variables
    int community = tm[u][f];
    int topic = zc[u][f];
    int m = ptrndata->users[u]->fol->members[f];
    nmc[m][community] -= 1;
    nrt[u][topic] -= 1;
    nmcsum[community] -= 1;
    nrtsum[u] -= 1;
    
    double Cgamma = C * gamma;
    double Kalpha = K * alpha;
    double Mdelta = M * delta;
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
        p_tc[k] = 0;
        for (int c = 0; c < C; c++) {
            p_tc[k] += (nmc[m][c] + delta) / (nmcsum[c] + Mdelta) *
            (nct[c][k] + gamma) / (nctsum[k] + Cgamma) *
            (ndt[u][k] + nrt[u][k] + alpha) / (ndtsum[u] + nrtsum[u] + Kalpha);
        }
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
        p_tc[k] += p_tc[k - 1];
    }
    // scaled sample because of unnormalized p_tw[]
    double u_value = ((double)random() / RAND_MAX) * p_tc[K - 1];
    
    for (topic = 0; topic < K; topic++) {
        if (p_tc[topic] > u_value) {
            break;
        }
    }
    
    // add newly estimated zc_i to count variables
    nmc[m][community] += 1;
    nrt[u][topic] += 1;
    nmcsum[community] += 1;
    nrtsum[u] += 1;
    return topic;
}

int model::sampling_tm(int u, int f) {
    // remove tm_i from the count variables
    int community = tm[u][f];
    int m = ptrndata->users[u]->fol->members[f];
    nmc[m][community] -= 1;
    nrc[u][community] -= 1;
    nmcsum[community] -= 1;
    nrcsum[u] -= 1;
    
    double Cgamma = C * gamma;
    double Mdelta = M * delta;
    double Kalpha = K * alpha;
    // do multinomial sampling via cumulative method
    for (int c = 0; c < C; c++) {
        for (int k = 0; k < K; k++) {
            p_cm[c] = (nmc[m][c] + delta) / (nmcsum[c] + Mdelta) *
            (nct[c][k] + gamma) / (nctsum[k] + Cgamma)*
            (ndt[u][k] + nrt[u][k] + alpha) / (ndtsum[u] + nrtsum[u] + Kalpha);
        }
    }
    // cumulate multinomial parameters
    for (int c = 1; c < C; c++) {
        p_cm[c] += p_cm[c - 1];
    }
    // scaled sample because of unnormalized p_tw[]
    double u_value = ((double)random() / RAND_MAX) * p_cm[C - 1];
    
    for (community = 0; community < C; community++) {
        if (p_cm[community] > u_value) {
            break;
        }
    }
    
    // add newly estimated zw_i to count variables
    nmc[m][community] += 1;
    nrc[u][community] += 1;
    nmcsum[community] += 1;
    nrcsum[u] += 1;
    
    return community;
}

void model::compute_theta() {
    for (int u = 0; u < U; u++) {
        for (int k = 0; k < K; k++) {
            theta[u][k] = (ndt[u][k] + alpha) / (ndtsum[u] + K * alpha);
        }
    }
}

void model::compute_phi() {
    for (int k = 0; k < K; k++) {
        for (int w = 0; w < V; w++) {
            phi[k][w] = (nwt[w][k] + beta) / (nwtsum[k] + V * beta);
        }
    }
}

void model::compute_lambda() {
    for (int k = 0; k < K; k++) {
        for (int c = 0; c < C; c++) {
            lambda[k][c] = (nct[c][k] + gamma) / (nctsum[k] + C * gamma);
        }
    }
}

void model::compute_omega() {
    for (int c = 0; c < C; c++) {
        for (int m = 0; m < M; m++) {
            omega[c][m] = (nmc[m][c] + delta) / (nmcsum[c] + M * delta);
        }
    }
}
/*
int model::load_model(string model_name) {
    int i, j;
    
    string filename = dir + model_name + tassignword_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %d to load model!\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_LONG];
    string line;

    // allocate memory for z and ptrndata
    z = new int*[M];
    ptrndata = new dataset(M);
    ptrndata->V = V;

    for (i = 0; i < M; i++) {
	char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
	if (!pointer) {
	    printf("Invalid word-topic assignment file, check the number of docs!\n");
	    return 1;
	}
	
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> words;
	vector<int> topics;
	for (j = 0; j < length; j++) {
	    string token = strtok.token(j);
    
	    strtokenizer tok(token, ":");
	    if (tok.count_tokens() != 2) {
		printf("Invalid word-topic assignment line!\n");
		return 1;
	    }
	    
	    words.push_back(atoi(tok.token(0).c_str()));
	    topics.push_back(atoi(tok.token(1).c_str()));
	}
	
	// allocate and add new document to the corpus
	document * pdoc = new document(words);
	ptrndata->add_doc(pdoc, i);
	
	// assign values for z
	z[i] = new int[topics.size()];
	for (j = 0; j < topics.size(); j++) {
	    z[i][j] = topics[j];
	}
    }   
    
    fclose(fin);
    
    return 0;
}
*/
int model::save_model(string model_name) {
    if (save_model_tassignword(dir + model_name + tassignword_suffix)) {
	return 1;
    }
    
    if (save_model_tassigncommunity(dir + model_name + tassigncommunity_suffix)) {
        return 1;
    }
    
    if (save_model_cassignmember(dir + model_name + cassignmember_suffix)) {
        return 1;
    }
    
    if (save_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_model_theta(dir + model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_model_phi(dir + model_name + phi_suffix)) {
	return 1;
    }
    
    if (save_model_lambda(dir + model_name + lambda_suffix)) {
        return 1;
    }
    
    if (save_model_omega(dir + model_name + omega_suffix)) {
        return 1;
    }
    
    if (twords > 0) {
	if (save_model_twords(dir + model_name + twords_suffix)) {
	    return 1;
	}
    }
    
    if (tcommunities > 0) {
        if (save_model_tcommunities(dir + model_name + tcommunities_suffix))
        {
            return 1;
        }
    }
    
    if (cmembers > 0) {
        if (save_model_cmembers(dir + model_name + cmembers_suffix))
        {
            return 1;
        }
    }
    
    return 0;
}

int model::save_model_tassignword(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->U; i++) {
	for (j = 0; j < ptrndata->users[i]->doc->length; j++) {
	    fprintf(fout, "%d:%d ", ptrndata->users[i]->doc->words[j], zw[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_model_tassigncommunity(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    
    // wirte docs with topic assignments for communities
    for (i = 0; i < ptrndata->U; i++) {
        for (j = 0; j < ptrndata->users[i]->fol->length; j++) {
            fprintf(fout, "%d:%d ", ptrndata->users[i]->fol->members[j], zc[i][j]);
        }
        fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_cassignmember(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    
    // wirte docs with topic assignments for communities
    for (i = 0; i < ptrndata->U; i++) {
        for (j = 0; j < ptrndata->users[i]->fol->length; j++) {
            fprintf(fout, "%d:%d ", ptrndata->users[i]->fol->members[j], tm[i][j]);
        }
        fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_theta(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < U; i++) {
	for (int j = 0; j < K; j++) {
	    fprintf(fout, "%f ", theta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_phi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < V; j++) {
	    fprintf(fout, "%f ", phi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_lambda(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < C; j++) {
            fprintf(fout, "%f ", lambda[i][j]);
        }
        fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_omega(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    
    for (int i = 0; i < C; i++) {
        for (int j = 0; j < M; j++) {
            fprintf(fout, "%f ", omega[i][j]);
        }
        fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "gamma=%f\n", gamma);
    fprintf(fout, "delta=%f\n", delta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "nusers=%d\n", U);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "nmembers=%d\n", M);
    fprintf(fout, "liter=%d\n", liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > V) {
	twords = V;
    }
    mapid2word::iterator it;
    
    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < V; w++) {
	    word_prob.first = w;
	    word_prob.second = phi[k][w];
	    words_probs.push_back(word_prob);
	}
    
        // quick sort to sort word-topic probability
	utils::quicksort(words_probs, 0, words_probs.size() - 1);
	
	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    it = id2word.find(words_probs[i].first);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }
    
    fclose(fout);    
    
    return 0;    
}

int model::save_model_tcommunities(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    
    if (tcommunities > C) {
        tcommunities = C;
    }
    
    for (int k = 0; k < K; k++) {
        vector<pair<int, double> > communities_probs;
        pair<int, double> community_prob;
        for (int c = 0; c < C; c++) {
            community_prob.first = c;
            community_prob.second = lambda[k][c];
            communities_probs.push_back(community_prob);
        }
        
        // quick sort to sort word-topic probability
        utils::quicksort(communities_probs, 0, communities_probs.size() - 1);
        
        fprintf(fout, "Topic %dth:\n", k);
        for (int i = 0; i < tcommunities; i++) {
            fprintf(fout, "\t%d   %f\n", communities_probs[i].first, communities_probs[i].second);
        }
    }
    
    fclose(fout);    
    
    return 0;    
}

int model::save_model_cmembers(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to save!\n", filename.c_str());
        return 1;
    }
    
    if (cmembers > M) {
        cmembers = M;
    }
    mapid2member::iterator it;
    
    for (int c = 0; c < C; c++) {
        vector<pair<int, double> > members_probs;
        pair<int, double> member_prob;
        for (int m = 0; m < M; m++) {
            member_prob.first = m;
            member_prob.second = omega[c][m];
            members_probs.push_back(member_prob);
        }
        
        // quick sort to sort word-topic probability
        utils::quicksort(members_probs, 0, members_probs.size() - 1);
        
        fprintf(fout, "Community %dth:\n", c);
        for (int i = 0; i < cmembers; i++) {
            it = id2member.find(members_probs[i].first);
            if (it != id2member.end()) {
                fprintf(fout, "\t%s   %f\n", (it->second).c_str(), members_probs[i].second);
            }
        }
    }
    
    fclose(fout);    
    
    return 0;    
}

/*
void model::inference() {
    if (twords > 0) {
        // print out top words per topic
        dataset::read_wordmap(dir + wordmapfile, &id2word);
    }
    
    printf("Sampling %d iterations for inference!\n", niters);
    
    for (inf_liter = 1; inf_liter <= niters; inf_liter++) {
        printf("Iteration %d ...\n", inf_liter);
        
        // for all newz_i
        for (int m = 0; m < newM; m++) {
            for (int n = 0; n < pnewdata->docs[m]->length; n++) {
                // (newz_i = newz[m][n])
                // sample from p(z_i|z_-i, w)
                int topic = inf_sampling(m, n);
                newz[m][n] = topic;
            }
        }
    }
    
    printf("Gibbs sampling for inference completed!\n");
    printf("Saving the inference outputs!\n");
    compute_newtheta();
    compute_newphi();
    inf_liter--;
    save_inf_model(cfile);
}

int model::inf_sampling_zw(int m, int n) {
    // remove z_i from the count variables
    int topic = newz[m][n];
    int w = pnewdata->docs[m]->words[n];
    int _w = pnewdata->_docs[m]->words[n];
    newnw[_w][topic] -= 1;
    newnd[m][topic] -= 1;
    newnwsum[topic] -= 1;
    newndsum[m] -= 1;
    
    double Vbeta = V * beta;
    double Kalpha = K * alpha;
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
        p[k] = (nw[w][k] + newnw[_w][k] + beta) / (nwsum[k] + newnwsum[k] + Vbeta) *
        (newnd[m][k] + alpha) / (newndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
        p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
        if (p[topic] > u) {
            break;
        }
    }
    
    // add newly estimated z_i to count variables
    newnw[_w][topic] += 1;
    newnd[m][topic] += 1;
    newnwsum[topic] += 1;
    newndsum[m] += 1;
    
    return topic;
}

void model::compute_newtheta() {
    for (int m = 0; m < newM; m++) {
        for (int k = 0; k < K; k++) {
            newtheta[m][k] = (newnd[m][k] + alpha) / (newndsum[m] + K * alpha);
        }
    }
}

void model::compute_newphi() {
    map<int, int>::iterator it;
    for (int k = 0; k < K; k++) {
        for (int w = 0; w < newV; w++) {
            it = pnewdata->_id2id.find(w);
            if (it != pnewdata->_id2id.end()) {
                newphi[k][w] = (nw[it->second][k] + newnw[w][k] + beta) / (nwsum[k] + newnwsum[k] + V * beta);
            }
        }
    }
}


int model::save_inf_model(string model_name) {
    if (save_inf_model_tassignword(dir + model_name + tassignword_suffix)) {
	return 1;
    }
    
    if (save_inf_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newphi(dir + model_name + phi_suffix)) {
	return 1;
    }

    if (twords > 0) {
	if (save_inf_model_twords(dir + model_name + twords_suffix)) {
	    return 1;
	}
    }
    
    return 0;
}

int model::save_inf_model_tassignword(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < pnewdata->M; i++) {    
	for (j = 0; j < pnewdata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", pnewdata->docs[i]->words[j], newz[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newtheta(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (i = 0; i < newM; i++) {
	for (j = 0; j < K; j++) {
	    fprintf(fout, "%f ", newtheta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newphi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < newV; j++) {
	    fprintf(fout, "%f ", newphi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "nusers=%d\n", newU);
    fprintf(fout, "nwords=%d\n", newV);
    fprintf(fout, "liter=%d\n", inf_liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > newV) {
	twords = newV;
    }
    mapid2word::iterator it;
    map<int, int>::iterator _it;
    
    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < newV; w++) {
	    word_prob.first = w;
	    word_prob.second = newphi[k][w];
	    words_probs.push_back(word_prob);
	}
    
        // quick sort to sort word-topic probability
	utils::quicksort(words_probs, 0, words_probs.size() - 1);
	
	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    _it = pnewdata->_id2id.find(words_probs[i].first);
	    if (_it == pnewdata->_id2id.end()) {
		continue;
	    }
	    it = id2word.find(_it->second);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }
    
    fclose(fout);    
    
    return 0;    
}
*/
