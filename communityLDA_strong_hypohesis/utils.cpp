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

#include <stdio.h>
#include <string>
#include <map>
#include "strtokenizer.h"
#include "utils.h"
#include "model.h"

using namespace std;

int utils::parse_args(int argc, char ** argv, model * pmodel) {
    int model_status = MODEL_STATUS_UNKNOWN;
    string dir = "";
    string model_name = "";
    string cfile = "";
    string ffile = "";
    double alpha = -1.0;
    double beta = -1.0;
    double gamma = -1.0;
    double delta = -1.0;
    int K = 0;
    int C = 0;
    int niters = 0;
    int savestep = 0;
    int twords = 0;
    int tcommunities = 0;
    int cmembers = 0;
    int withrawdata = 0;

    int i = 0; 
    while (i < argc) {
	string arg = argv[i];
	
	if (arg == "-est") {
	    model_status = MODEL_STATUS_EST;
	    
	} else if (arg == "-estc") {
	    model_status = MODEL_STATUS_ESTC;
	    
	} else if (arg == "-inf") {
	    model_status = MODEL_STATUS_INF;
	    
	} else if (arg == "-dir") {
	    dir = argv[++i];	    
	    
	} else if (arg == "-cfile") {
	    cfile = argv[++i];	    
	    
	} else if (arg == "-ffile") {
        ffile = argv[++i];
    } else if (arg == "-model") {
	    model_name = argv[++i];	    	    
	    
	} else if (arg == "-alpha") {
	    alpha = atof(argv[++i]);	    
	    
	} else if (arg == "-beta") {
	    beta = atof(argv[++i]);	    
	    
	} else if (arg == "-gamma") {
	    gamma = atof(argv[++i]);
	    
	} else if (arg == "-delta") {
	    delta = atof(argv[++i]);
	    
	} else if (arg == "-ntopics") {
	    K = atoi(argv[++i]);	    
	    
	} else if (arg == "-ncommunities") {
	    C = atoi(argv[++i]);
	    
	} else if (arg == "-niters") {
	    niters = atoi(argv[++i]);	    
	    
	} else if (arg == "-savestep") {
	    savestep = atoi(argv[++i]);
	    
	} else if (arg == "-twords") {
	    twords = atoi(argv[++i]);
	    
	} else if (arg == "-tcommunities") {
	    tcommunities = atoi(argv[++i]);
	    
	} else if (arg == "-cmembers") {
	    cmembers = atoi(argv[++i]);
	    
	} else if (arg == "-withrawdata") {
	    withrawdata = 1;
	
	} else {
	    // any more?
	}	
		
	i++;
    }
    
    if (model_status == MODEL_STATUS_EST) {
        if (cfile == "" || ffile == "") {
            printf("Please specify the input data file for model estimation!\n");
            return 1;
        }
	
        pmodel->model_status = model_status;
	
        if (K > 0) {
            pmodel->K = K;
        }
    
        if (C > 0) {
            pmodel->C = C;
        }
	
        if (alpha >= 0.0) {
            pmodel->alpha = alpha;
        } else {
            // default value for alpha
            pmodel->alpha = 50.0 / pmodel->K;
        }
	
        if (beta >= 0.0) {
            pmodel->beta = beta;
        }
        
        if (gamma >= 0.0) {
            pmodel->gamma = gamma;
        }
        
        if (delta >= 0.0) {
            pmodel->delta = delta;
        }
	
        if (niters > 0) {
            pmodel->niters = niters;
        }
	
        if (savestep > 0) {
            pmodel->savestep = savestep;
        }
	
        if (twords > 0) {
            pmodel->twords = twords;
        }
        
        if (tcommunities > 0) {
            pmodel->tcommunities = tcommunities;
        }
        
        if (cmembers > 0) {
            pmodel->cmembers = cmembers;
        }
	
        pmodel->cfile = cfile;
        pmodel->ffile = ffile;
	
        //the content file and the follow file need be in the same dir
        string::size_type idx = ffile.find_last_of("/");
        if (idx == string::npos) { //if not find
            dir = "./";
        } else {
            dir = ffile.substr(0, idx + 1);
        }
        idx = cfile.find_last_of("/");
        if (idx == string::npos) { //if not find
            pmodel->dir = "./";
            if (pmodel->dir != dir) {
                printf("the content file and the follow file need be in the same dir");
                return 1;
            }
        } else {
            pmodel->dir = cfile.substr(0, idx + 1);
            if (pmodel->dir != dir) {
                printf("the content file and the follow file need be in the same dir");
                return 1;
            }
            pmodel->cfile = cfile.substr(idx + 1, cfile.size() - pmodel->dir.size());
            pmodel->ffile = ffile.substr(idx + 1, ffile.size() - pmodel->dir.size());
            printf("dir = %s\n", pmodel->dir.c_str());
            printf("cfile = %s\n", pmodel->cfile.c_str());
            printf("ffile = %s\n", pmodel->ffile.c_str());
        }
    } 
    
    if (model_status == MODEL_STATUS_ESTC) {
        if (dir == "") {
            printf("Please specify model directory!\n");
            return 1;
        }
	
        if (model_name == "") {
            printf("Please specify model name upon that you want to continue estimating!\n");
            return 1;
        }	

        pmodel->model_status = model_status;

        if (dir[dir.size() - 1] != '/') {
            dir += "/";
        }
        pmodel->dir = dir;

        pmodel->model_name = model_name;

        if (niters > 0) {
            pmodel->niters = niters;
        }
	
        if (savestep > 0) {
            pmodel->savestep = savestep;
        }
	
        if (twords > 0) {
            pmodel->twords = twords;
        }
        
        if (tcommunities > 0) {
            pmodel->tcommunities = tcommunities;
        }
        
        if (cmembers > 0) {
            pmodel->cmembers = cmembers;
        }
	
        // read <model>.others file to assign values for ntopics, alpha, beta, etc.
        if (read_and_parse(pmodel->dir + pmodel->model_name + pmodel->others_suffix, pmodel)) {
            return 1;
        }	
    } 
    
    if (model_status == MODEL_STATUS_INF) {
        if (dir == "") {
            printf("Please specify model directory please!\n");
            return 1;
        }
	
        if (model_name == "") {
            printf("Please specify model name for inference!\n");
            return 1;
        }	

        if (cfile == "" || ffile == "") {
            printf("Please specify the new data file for inference!\n");
            return 1;
        }
	
        pmodel->model_status = model_status;

        if (dir[dir.size() - 1] != '/') {
            dir += "/";
        }
        pmodel->dir = dir;
	
        pmodel->model_name = model_name;

        pmodel->cfile = cfile;
        pmodel->ffile = ffile;

        if (niters > 0) {
            pmodel->niters = niters;
        } else {
            // default number of Gibbs sampling iterations for doing inference
            pmodel->niters = 20;
        }
	
        if (twords > 0) {
            pmodel->twords = twords;
        }
        
        if (tcommunities > 0) {
            pmodel->tcommunities = tcommunities;
        }
        
        if (cmembers > 0) {
            pmodel->cmembers = cmembers;
        }
	
        if (withrawdata > 0) {
            pmodel->withrawstrs = withrawdata;
        }
		
        // read <model>.others file to assign values for ntopics, alpha, beta, etc.
        if (read_and_parse(pmodel->dir + pmodel->model_name + pmodel->others_suffix, pmodel)) {
            return 1;
        }
    }
    
    if (model_status == MODEL_STATUS_UNKNOWN) {
	printf("Please specify the task you would like to perform (-est/-estc/-inf)!\n");
	return 1;
    }
    
    return 0;
}

int utils::read_and_parse(string filename, model * pmodel) {
    // open file <model>.others to read:
    // alpha=?
    // beta=?
    // gamma=?
    // delta=?
    // ntopics=?
    // ncommunities=?
    // nusers=?
    // nwords=?
    // nmemebers=?
    // liter=? // current iteration (when the model was saved)
    
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file: %s\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    while (fgets(buff, BUFF_SIZE_SHORT - 1, fin)) {
	line = buff;
	strtokenizer strtok(line, "= \t\r\n");
	int count = strtok.count_tokens();
	
	if (count != 2) {
	    // invalid, ignore this line
	    continue;
	}

	string optstr = strtok.token(0);
	string optval = strtok.token(1);
	
	if (optstr == "alpha") {
	    pmodel->alpha = atof(optval.c_str());
	    
	} else if (optstr == "beta") {	    
	    pmodel->beta = atof(optval.c_str());
	
	} else if (optstr == "gamma") {
	    pmodel->gamma = atof(optval.c_str());
        
	} else if (optstr == "delta") {
	    pmodel->delta = atof(optval.c_str());
        
	} else if (optstr == "ntopics") {
	    pmodel->K = atoi(optval.c_str());
	
	} else if (optstr == "ncommunities") {
	    pmodel->C = atoi(optval.c_str());
        
	} else if (optstr == "nusers") {
	    pmodel->U = atoi(optval.c_str());
	 
	} else if (optstr == "nwords") {
	    pmodel->V = atoi(optval.c_str());
	
	} else if (optstr == "nmembers") {
	    pmodel->M = atoi(optval.c_str());
        
	} else if (optstr == "liter") {
	    pmodel->liter = atoi(optval.c_str());
	
	} else {
	    // any more?
	}
    }
    
    fclose(fin);
    
    return 0;
}

string utils::generate_model_name(int iter) {
    string model_name = "model-";

    char buff[BUFF_SIZE_SHORT];
    
    if (0 <= iter && iter < 10) {
	sprintf(buff, "0000%d", iter);
    } else if (10 <= iter && iter < 100) {
	sprintf(buff, "000%d", iter);
    } else if (100 <= iter && iter < 1000) {
	sprintf(buff, "00%d", iter);
    } else if (1000 <= iter && iter < 10000) {
	sprintf(buff, "0%d", iter);
    } else {
	sprintf(buff, "%d", iter);
    }
    
    if (iter >= 0) {
	model_name += buff;
    } else {
	model_name += "final";
    }
    
    return model_name;
}

void utils::sort(vector<double> & probs, vector<int> & words) {
    for (int i = 0; i < probs.size() - 1; i++) {
	for (int j = i + 1; j < probs.size(); j++) {
	    if (probs[i] < probs[j]) {
		double tempprob = probs[i];
		int tempword = words[i];
		probs[i] = probs[j];
		words[i] = words[j];
		probs[j] = tempprob;
		words[j] = tempword;
	    }
	}
    }
}

void utils::quicksort(vector<pair<int, double> > & vect, int left, int right) {
    int l_hold, r_hold;
    pair<int, double> pivot;
    
    l_hold = left;
    r_hold = right;    
    int pivotidx = left;
    pivot = vect[pivotidx];

    while (left < right) {
	while (vect[right].second <= pivot.second && left < right) {
	    right--;
	}
	if (left != right) {
	    vect[left] = vect[right];
	    left++;
	}
	while (vect[left].second >= pivot.second && left < right) {
	    left++;
	}
	if (left != right) {
	    vect[right] = vect[left];
	    right--;
	}
    }

    vect[left] = pivot;
    pivotidx = left;
    left = l_hold;
    right = r_hold;
    
    if (left < pivotidx) {
	quicksort(vect, left, pivotidx - 1);
    }
    if (right > pivotidx) {
	quicksort(vect, pivotidx + 1, right);
    }    
}

