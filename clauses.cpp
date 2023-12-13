/*
Dynex SDK (beta) Neuromorphic Computing Library
Copyright (c) 2021-2023, Dynex Developers

All rights reserved.

1. Redistributions of source code must retain the above copyright notice, this list of
    conditions and the following disclaimer.
 
2. Redistributions in binary form must reproduce the above copyright notice, this list
   of conditions and the following disclaimer in the documentation and/or other
   materials provided with the distribution.
 
3. Neither the name of the copyright holder nor the names of its contributors may be
   used to endorse or promote products derived from this software without specific
   prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "clauses.hpp"

ULL MAXWEIGHT = (1ULL<<63) - 1;

//! doubles the amount of memory available to store clauses
void Clauses::doubleStorage() {
	capacity <<= 1;
	++logcapacity;
	int_c *clauses2 = new int_c[capacity];
	assert(clauses2 != NULL);
	memcpy(clauses2, clauses, sizeof(int_c) * capacity/2);
	delete [] clauses;
	clauses = clauses2;
	int *headsFreeList2 = new int[logcapacity];
	assert(headsFreeList2 != NULL);
	memcpy(headsFreeList2, headsFreeList, sizeof(int) * (logcapacity-1));
	delete [] headsFreeList;
	headsFreeList = headsFreeList2;
	assert(headsFreeList[logcapacity-2] == 0);
	// add the new storage block to the buddy memory management
	headsFreeList[logcapacity-2] = capacity/2;
	headsFreeList[logcapacity-1] = 0;
	clauses[capacity/2] = 0;
}

//! Clauses constructor
Clauses::Clauses() {
	assigned = NULL;
	logcapacity = 16;
	capacity = 1<<logcapacity;
	clauses = new int_c[capacity];
	headsFreeList = new int[logcapacity];
	for (int i=0; i+1<logcapacity; ++i) {
		headsFreeList[i] = 1<<i;
		clauses[1<<i] = 0;
	}
	headsFreeList[logcapacity-1] = 0;
}

//! Clauses destructor
Clauses::~Clauses() {
	delete [] clauses;
	assert(assigned != NULL);
	assigned -= nVars;
	delete [] assigned;
	delete [] headsFreeList;
}
