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
#include <stdlib.h>
#include "constants.h"
#include "strtokenizer.h"
#include "dataset.h"

using namespace std;

int dataset::write_wordmap(string wordmapfile, mapword2id * pword2id) {
    FILE * fout = fopen(wordmapfile.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to write!\n", wordmapfile.c_str());
	return 1;
    }    
    
    mapword2id::iterator it;
    fprintf(fout, "%d\n", pword2id->size());
    for (it = pword2id->begin(); it != pword2id->end(); it++) {
	fprintf(fout, "%s %d\n", (it->first).c_str(), it->second);
    }
    
    fclose(fout);
    
    return 0;
}

int dataset::write_membermap(string membermapfile, mapmember2id * pmember2id) {
    FILE * fout = fopen(membermapfile.c_str(), "w");
    if (!fout) {
        printf("Cannot open file %s to write!\n", membermapfile.c_str());
        return 1;
    }
    
    mapmember2id::iterator it;
    fprintf(fout, "%d\n", pmember2id->size());
    for (it = pmember2id->begin(); it != pmember2id->end(); it++) {
        fprintf(fout, "%s %d\n", (it->first).c_str(), it->second);
    }
    
    fclose(fout);
    
    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapword2id * pword2id) {
    pword2id->clear();
    
    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);
    
    for (int i = 0; i < nwords; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, " \t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}
	
	pword2id->insert(pair<string, int>(strtok.token(0), atoi(strtok.token(1).c_str())));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapid2word * pid2word) {
    pid2word->clear();
    
    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);
    
    for (int i = 0; i < nwords; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, " \t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}
	
	pid2word->insert(pair<int, string>(atoi(strtok.token(1).c_str()), strtok.token(0)));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_membermap(string membermapfile, mapmember2id * pmember2id) {
    pmember2id->clear();
    
    FILE * fin = fopen(membermapfile.c_str(), "r");
    if (!fin) {
        printf("Cannot open file %s to read!\n", membermapfile.c_str());
        return 1;
    }
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nmembers = atoi(buff);
    
    for (int i = 0; i < nmembers; i++) {
        fgets(buff, BUFF_SIZE_SHORT - 1, fin);
        line = buff;
        
        strtokenizer strtok(line, " \t\r\n");
        if (strtok.count_tokens() != 2) {
            continue;
        }
        
        pmember2id->insert(pair<string, int>(strtok.token(0), atoi(strtok.token(1).c_str())));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_membermap(string membermapfile, mapid2member * pid2member) {
    pid2member->clear();
    
    FILE * fin = fopen(membermapfile.c_str(), "r");
    if (!fin) {
        printf("Cannot open file %s to read!\n", membermapfile.c_str());
        return 1;
    }
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nmembers = atoi(buff);
    
    for (int i = 0; i < nmembers; i++) {
        fgets(buff, BUFF_SIZE_SHORT - 1, fin);
        line = buff;
        
        strtokenizer strtok(line, " \t\r\n");
        if (strtok.count_tokens() != 2) {
            continue;
        }
        
        pid2member->insert(pair<int, string>(atoi(strtok.token(1).c_str()), strtok.token(0)));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_trndata(string cfile, string ffile, string wordmapfile, string membermapfile) {
    mapword2id word2id;
    mapmember2id member2id;
    
    FILE * cfin = fopen(cfile.c_str(), "r");
    FILE * ffin = fopen(ffile.c_str(), "r");
    if (!cfin) {
	printf("Cannot open file %s to read!\n", cfile.c_str());
	return 1;
    }
    if (!ffin) {
        printf("Cannot open file %s to read!\n", ffile.c_str());
        return 1;
    }
    
    mapword2id::iterator wit;
    mapmember2id::iterator mit;
    
    char cbuff[BUFF_SIZE_LONG];
    char fbuff[BUFF_SIZE_LONG];
    string cline;
    string fline;
    
    // get the number of users
    fgets(cbuff, BUFF_SIZE_LONG - 1, cfin);
    fgets(fbuff, BUFF_SIZE_LONG - 1, ffin);
    if ((U = atoi(cbuff)) != atoi(fbuff)) {
        printf("the number of users in contentfile and followfile unequal!");
        return 1;
    }

    if (U <= 0) {
	printf("No document available!\n");
	return 1;
    }
    
    // allocate memory for corpus
    if (users) {
	deallocate();
    } else {
	users = new user*[U];
    }
    
    // set number of words to zero
    V = 0;
    M = 0;
    
    for (int i = 0; i < U; i++) {
        fgets(cbuff, BUFF_SIZE_LONG - 1, cfin);
        fgets(fbuff, BUFF_SIZE_LONG - 1, ffin);
        cline = cbuff;
        fline = fbuff;
        strtokenizer cstrtok(cline, " \t\r\n");
        int clength = cstrtok.count_tokens();
        strtokenizer fstrtok(fline, " \t\r\n");
        int flength = fstrtok.count_tokens();
        if (clength <= 0 || flength <=0) {
            printf("Invalid (empty) document!\n");
            deallocate();
            U = V = M = 0;
            return 1;
        }
	
        // allocate new user
        user * puser = new user(clength, flength);
	
        for (int j = 0; j < clength; j++) {
            wit = word2id.find(cstrtok.token(j));
            if (wit == word2id.end()) {
                // word not found, i.e., new word
                puser->doc->words[j] = word2id.size();
                word2id.insert(pair<string, int>(cstrtok.token(j), word2id.size()));
            } else {
                puser->doc->words[j] = wit->second;
            }
        }
        
        for (int k = 0; k < flength; k++) {
            mit = member2id.find(fstrtok.token(k));
            if (mit == member2id.end()) {
                // word not found, i.e., new word
                puser->fol->members[k] = member2id.size();
                member2id.insert(pair<string, int>(fstrtok.token(k), member2id.size()));
            } else {
                puser->fol->members[k] = mit->second;
            }
        }
	
        // add new doc to the corpus
        add_user(puser, i);
    }
    
    fclose(cfin);
    fclose(ffin);
    
    // write word map to file
    if (write_wordmap(wordmapfile, &word2id)) {
	return 1;
    }
    
    if (write_membermap(membermapfile, &member2id)) {
        return 1;
    }
    
    // update number of words
    V = word2id.size();
    M = member2id.size();
    
    return 0;
}
/*
int dataset::read_newdata(string cfile, string wordmapfile) {
    mapword2id word2id;
    map<int, int> id2_id;
    
    read_wordmap(wordmapfile, &word2id);
    if (word2id.size() <= 0) {
	printf("No word map available!\n");
	return 1;
    }

    FILE * fin = fopen(cfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", cfile.c_str());
	return 1;
    }   

    mapword2id::iterator it;
    map<int, int>::iterator _it;
    char buff[BUFF_SIZE_LONG];
    string line;
    
    // get number of new documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	printf("No document available!\n");
	return 1;
    }
    
    // allocate memory for corpus
    if (users) {
	deallocate();
    } else {
	users = new user*[U];
    }
    _users = new user*[M];
    
    // set number of words to zero
    V = 0;
    
    for (int i = 0; i < M; i++) {
	fgets(buff, BUFF_SIZE_LONG - 1, fin);
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> doc;
	vector<int> _doc;
	for (int j = 0; j < length; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., word unseen in training data
		// do anything? (future decision)
	    } else {
		int _id;
		_it = id2_id.find(it->second);
		if (_it == id2_id.end()) {
		    _id = id2_id.size();
		    id2_id.insert(pair<int, int>(it->second, _id));
		    _id2id.insert(pair<int, int>(_id, it->second));
		} else {
		    _id = _it->second;
		}
		
		doc.push_back(it->second);
		_doc.push_back(_id);
	    }
	}
	
	// allocate memory for new doc
	document * pdoc = new document(doc);
	document * _pdoc = new document(_doc);
	
	// add new doc
	add_doc(pdoc, i);
	_add_doc(_pdoc, i);
    }
    
    fclose(fin);
    
    // update number of new words
    V = id2_id.size();
    
    return 0;
}
*/
/*
int dataset::read_newdata_withrawstrs(string cfile, string wordmapfile) {
    mapword2id word2id;
    map<int, int> id2_id;
    
    read_wordmap(wordmapfile, &word2id);
    if (word2id.size() <= 0) {
	printf("No word map available!\n");
	return 1;
    }

    FILE * fin = fopen(cfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", cfile.c_str());
	return 1;
    }   

    mapword2id::iterator it;
    map<int, int>::iterator _it;
    char buff[BUFF_SIZE_LONG];
    string line;
    
    // get number of new documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	printf("No document available!\n");
	return 1;
    }
    
    // allocate memory for corpus
    if (docs) {
	deallocate();
    } else {
	docs = new document*[M];
    }
    _docs = new document*[M];
    
    // set number of words to zero
    V = 0;
    
    for (int i = 0; i < M; i++) {
	fgets(buff, BUFF_SIZE_LONG - 1, fin);
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> doc;
	vector<int> _doc;
	for (int j = 0; j < length - 1; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., word unseen in training data
		// do anything? (future decision)
	    } else {
		int _id;
		_it = id2_id.find(it->second);
		if (_it == id2_id.end()) {
		    _id = id2_id.size();
		    id2_id.insert(pair<int, int>(it->second, _id));
		    _id2id.insert(pair<int, int>(_id, it->second));
		} else {
		    _id = _it->second;
		}
		
		doc.push_back(it->second);
		_doc.push_back(_id);
	    }
	}
	
	// allocate memory for new doc
	document * pdoc = new document(doc, line);
	document * _pdoc = new document(_doc, line);
	
	// add new doc
	add_doc(pdoc, i);
	_add_doc(_pdoc, i);
    }
    
    fclose(fin);
    
    // update number of new words
    V = id2_id.size();
    
    return 0;
}
*/
