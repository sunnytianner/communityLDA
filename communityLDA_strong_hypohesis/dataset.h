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

#ifndef	_DATASET_H
#define	_DATASET_H

#include <string>
#include <vector>
#include <map>

using namespace std;

// map of words/terms [string => int]
typedef map<string, int> mapword2id;
// map of words/terms [int => string]
typedef map<int, string> mapid2word;
// ##add## map of members [string => int]
typedef map<string, int> mapmember2id;
// ##add## map of members [int => string]
typedef map<int, string> mapid2member;

class document {
public:
    int * words;
    string rawstr;
    int length;
    
    document() {
	words = NULL;
	rawstr = "";
	length = 0;	
    }
    
    document(int length) {
	this->length = length;
	rawstr = "";
	words = new int[length];	
    }
    
    document(int length, int * words) {
	this->length = length;
	rawstr = "";
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = words[i];
	}
    }

    document(int length, int * words, string rawstr) {
	this->length = length;
	this->rawstr = rawstr;
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = words[i];
	}
    }
    
    document(vector<int> & doc) {
	this->length = doc.size();
	rawstr = "";
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = doc[i];
	}
    }

    document(vector<int> & doc, string rawstr) {
	this->length = doc.size();
	this->rawstr = rawstr;
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = doc[i];
	}
    }
    
    ~document() {
	if (words) {
	    delete words;
	}
    }
};

class follow {
public:
    int * members;
    int * communities;
    string rawstr;
    int length;
    
    follow() {
        members = NULL;
        communities = NULL;
        rawstr = "";
        length = 0;
    }
    
    follow(int length) {
        this->length = length;
        rawstr = "";
        members = new int[length];
        communities = new int[length];
    }
    
    follow(int length, int * members) {
        this->length = length;
        rawstr = "";
        this->members = new int[length];
        this->communities = new int[length];
        for (int i = 0; i < length; i++) {
            this->members[i] = members[i];
        }
    }
    
    follow(int length, int * members, string rawstr) {
        this->length = length;
        this->rawstr = rawstr;
        this->members = new int[length];
        this->communities = new int[length];
        for (int i = 0; i < length; i++) {
            this->members[i] = members[i];
        }
    }
    
    follow(vector<int> & fol) {
        this->length = fol.size();
        rawstr = "";
        this->members = new int[length];
        this->communities = new int[length];
        for (int i = 0; i < length; i++) {
            this->members[i] = fol[i];
        }
    }
    
    follow(vector<int> & fol, string rawstr) {
        this->length = fol.size();
        this->rawstr = rawstr;
        this->members = new int[length];
        this->communities = new int[length];
        for (int i = 0; i < length; i++) {
            this->members[i] = fol[i];
        }
    }
    
    ~follow() {
        if (members) {
            delete members;
        }
        if (communities) {
            delete communities;
        }
    }
};

class user {
public:
    document * doc;
    follow * fol;
    
    user() {
        doc = new document();
        fol = new follow();
    }
    
    user(int dlength, int flength) {
        doc = new document(dlength);
        fol = new follow(flength);
    }
    
    ~user() {
        if (doc) {
            delete doc;
        }
        if (fol) {
            delete fol;
        }
    }
};

class dataset {
public:
    user ** users;   //U users
    user ** _users; // used only for inference
    map<int, int> _id2id; // also used only for inference
    int U; // number of users
    int V; // number of words
    int M; // number of members
    
    dataset() {
	users = NULL;
	_users = NULL;
	U = 0;
	V = 0;
    M = 0;
    }
    
    dataset(int U) {
	this->U = U;
	this->V = 0;
    this->M = 0;
	users = new user*[U];
	_users = NULL;
    }   
    
    ~dataset() {
	if (users) {
	    for (int i = 0; i < U; i++) {
        users[i]->~user();
	    }
	}
	delete users;
	
	if (_users) {
	    for (int i = 0; i < U; i++) {
		_users[i]->~user();
	    }
	}
	delete _users;
    }
    
    void deallocate() {
	if (users) {
	    for (int i = 0; i < U; i++) {
		users[i]->~user();
	    }
	}
	delete users;
	users = NULL;

	if (_users) {
	    for (int i = 0; i < U; i++) {
		_users[i]->~user();
	    }
	}
	delete _users;
	_users = NULL;
    }
    
    void add_user(user * u, int idx) {
	if (0 <= idx && idx < U) {
	    users[idx] = u;
	}
    }   
    
    void _add_user(user * u, int idx) {
	if (0 <= idx && idx < U) {
	    _users[idx] = u;
	}
    }       

    static int write_wordmap(string wordmapfile, mapword2id * pword2id);
    static int write_membermap(string membermapfile, mapmember2id * pmember2id);
    static int read_wordmap(string wordmapfile, mapword2id * pword2id);
    static int read_wordmap(string wordmapfile, mapid2word * pid2word);
    static int read_membermap(string membermapfile, mapmember2id * pmember2id);
    static int read_membermap(string membermapfile, mapid2member * pid2member);
    
    int read_trndata(string cfile, string ffile, string wordmapfile, string membermapfile);
    //int read_newdata(string cfile, string wordmapfile);
    //int read_newdata_withrawstrs(string cfile, string wordmapfile);
};

#endif

