/*
 * HashTable.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: xgo
 */

#include "HashTableLongKey.h"

HashTableLongKey::HashTableLongKey(QueryDataset * qDataset) {
	this->hashKeyLength = Config::iKmer;
	this->dataSet = qDataset;
	this->vviHashTable = NULL;
	this->hashTableSize = 0;
	if (Config::bSpeicialVersion || Config::bMappingVersion) {
		this->iBitOfReadId = 35;
	} else {
		this->iBitOfReadId = 40;
	}
	this->iMask = (1L << this->iBitOfReadId) - 1;
	this->maxSingleHashCollision = 0;
	this->numberOfHashCollision = 0;
}

HashTableLongKey::~HashTableLongKey() {
	if (this->vviHashTable != NULL) {
		for (size_t i = 0; i < this->vviHashTable->size(); i++) {
			vector<INT64> * temp = this->vviHashTable->at(i);
			if (temp != NULL) {
				temp->clear();
				delete temp;
			}
		}

		this->vviHashTable->clear();
		delete this->vviHashTable;
		this->vviHashTable = NULL;
	}
}

bool HashTableLongKey::createHashTables() {

	if (this->dataSet == NULL) {
		cout << "no data set" << endl;
		return false;
	} else if (this->vviHashTable != NULL) {
		cout << "Hash Table already exists." << endl;
		return false;
	} else {
		CLOCKSTART
		MEMORYSTART
		UINT64 temp = Config::iHashTableFactor;
		if (Config::bMappingVersion) {
			UINT64 iTotalLength = 0;
			for (UINT64 i = 0; i < dataSet->getNumberOfReads(); i++) {
				iTotalLength += dataSet->vQueryReads->at(i)->getRawReadLength();
			}
			hashTableSize = temp * (iTotalLength / Config::iKmer);
		} else {
			hashTableSize = temp * this->dataSet->getNumberOfReads();
		}
		hashTableSize = getPrimeLargerThanNumber(hashTableSize);
		setHashTableSizeAndInitialize(hashTableSize);
		MEMORYSTOP
		CLOCKSTOP
		return true;
	}
}

bool compareSequence(const string & _seq1, INT64 & _start1, const string & _seq2, INT64 & _start2, INT64 _length) {
	for (INT64 i = 0; i < _length; i++) {
		if (_seq1.at(_start1 + i) != _seq2.at(_start2 + i)) {
			return false;
		}
	}
	return true;
}

vector<INT64> * HashTableLongKey::getListOfReads(const string & _sequence, INT64 _start, INT64 _length) {
	UINT64 index = -1;
	UINT64 currentCollision = 0;
	index = this->hashFunction(_sequence, _start, _length);
	while (!this->isEmptyAt(index)) {
		vector<INT64>* readList = this->vviHashTable->at(index);
		INT64 data = readList->at(0);
		INT64 keyReadID = data & this->iMask;
		INT64 keyPosition = data >> this->iBitOfReadId;
		if (compareSequence(_sequence, _start, (*this->dataSet->getReadFromID(keyReadID)->getSequence()), keyPosition,
				_length)) {
			break;
		}
		currentCollision++;
		if (currentCollision > this->maxSingleHashCollision)
			return NULL;
		index = (index == hashTableSize - 1) ? 0 : index + 1; 	// Increment the index
	}

	return this->vviHashTable->at(index);
}

UINT64 HashTableLongKey::getPrimeLargerThanNumber(UINT64 number) {
	// Pre-computed list of prime number.
	UINT64 array[452] = { 3359, 498857, 1114523, 1180043, 1245227, 1310759, 1376447, 1442087, 1507379, 1573667, 1638899,
			1704023, 1769627, 1835027, 1900667, 1966127, 2031839, 2228483, 2359559, 2490707, 2621447, 2752679, 2883767,
			3015527, 3145739, 3277283, 3408323, 3539267, 3670259, 3801143, 3932483, 4063559, 4456643, 4718699, 4980827,
			5243003, 5505239, 5767187, 6029603, 6291563, 6553979, 6816527, 7079159, 7340639, 7602359, 7864799, 8126747,
			8913119, 9437399, 9962207, 10485767, 11010383, 11534819, 12059123, 12583007, 13107923, 13631819, 14156543,
			14680067, 15204467, 15729647, 16253423, 17825999, 18874379, 19923227, 20971799, 22020227, 23069447,
			24117683, 25166423, 26214743, 27264047, 28312007, 29360147, 30410483, 31457627, 32505983, 35651783,
			37749983, 39845987, 41943347, 44040383, 46137887, 48234623, 50331707, 52429067, 54526019, 56623367,
			58720307, 60817763, 62915459, 65012279, 71303567, 75497999, 79691867, 83886983, 88080527, 92275307,
			96470447, 100663439, 104858387, 109052183, 113246699, 117440699, 121635467, 125829239, 130023683, 142606379,
			150994979, 159383759, 167772239, 176160779, 184549559, 192938003, 201327359, 209715719, 218104427,
			226493747, 234882239, 243269639, 251659139, 260047367, 285215507, 301989959, 318767927, 335544323,
			352321643, 369100463, 385876703, 402654059, 419432243, 436208447, 452986103, 469762067, 486539519,
			503316623, 520094747, 570425399, 603979919, 637534763, 671089283, 704643287, 738198347, 771752363,
			805307963, 838861103, 872415239, 905971007, 939525143, 973079279, 1006633283, 1040187419, 1140852767,
			1207960679, 1275069143, 1342177379, 1409288183, 1476395699, 1543504343, 1610613119, 1677721667, 1744830587,
			1811940419, 1879049087, 1946157419, 2013265967, 2080375127, 2281701827, 2415920939, 2550137039, 2684355383,
			2818572539, 2952791147, 3087008663, 3221226167, 3355444187, 3489661079, 3623878823, 3758096939, 3892314659,
			4026532187, 4160749883, 4563403379, 4831838783, 5100273923, 5368709219, 5637144743, 5905580687, 6174015503,
			6442452119, 6710886467, 6979322123, 7247758307, 7516193123, 7784629079, 8053065599, 8321499203, 9126806147,
			9663676523, 10200548819, 10737418883, 11274289319, 11811160139, 12348031523, 12884902223, 13421772839,
			13958645543, 14495515943, 15032386163, 15569257247, 16106127887, 16642998803, 18253612127, 19327353083,
			20401094843, 21474837719, 22548578579, 23622320927, 24696062387, 25769803799, 26843546243, 27917287907,
			28991030759, 30064772327, 31138513067, 32212254947, 33285996803, 36507222923, 38654706323, 40802189423,
			42949673423, 45097157927, 47244640319, 49392124247, 51539607599, 53687092307, 55834576979, 57982058579,
			60129542339, 62277026327, 64424509847, 66571993199, 73014444299, 77309412407, 81604379243, 85899346727,
			90194314103, 94489281203, 98784255863, 103079215439, 107374183703, 111669150239, 115964117999, 120259085183,
			124554051983, 128849019059, 133143986399, 146028888179, 154618823603, 163208757527, 171798693719,
			180388628579, 188978561207, 197568495647, 206158430447, 214748365067, 223338303719, 231928234787,
			240518168603, 249108103547, 257698038539, 266287975727, 292057776239, 309237645803, 326417515547,
			343597385507, 360777253763, 377957124803, 395136991499, 412316861267, 429496730879, 446676599987,
			463856468987, 481036337207, 498216206387, 515396078039, 532575944723, 584115552323, 618475290887,
			652835029643, 687194768879, 721554506879, 755914244627, 790273985219, 824633721383, 858993459587,
			893353198763, 927712936643, 962072674643, 996432414899, 1030792152539, 1065151889507, 1168231105859,
			1236950582039, 1305670059983, 1374389535587, 1443109012607, 1511828491883, 1580547965639, 1649267441747,
			1717986918839, 1786706397767, 1855425872459, 1924145348627, 1992864827099, 2061584304323, 2130303780503,
			2336462210183, 2473901164367, 2611340118887, 2748779070239, 2886218024939, 3023656976507, 3161095931639,
			3298534883999, 3435973836983, 3573412791647, 3710851743923, 3848290698467, 3985729653707, 4123168604483,
			4260607557707, 4672924419707, 4947802331663, 5222680234139, 5497558138979, 5772436047947, 6047313952943,
			6322191860339, 6597069767699, 6871947674003, 7146825580703, 7421703488567, 7696581395627, 7971459304163,
			8246337210659, 8521215117407, 9345848837267, 9895604651243, 10445360463947, 10995116279639, 11544872100683,
			12094627906847, 12644383722779, 13194139536659, 13743895350023, 14293651161443, 14843406975659,
			15393162789503, 15942918604343, 16492674420863, 17042430234443, 18691697672867, 19791209300867,
			20890720927823, 21990232555703, 23089744183799, 24189255814847, 25288767440099, 26388279068903,
			27487790694887, 28587302323787, 29686813951463, 30786325577867, 31885837205567, 32985348833687,
			34084860462083, 37383395344739, 39582418600883, 41781441856823, 43980465111383, 46179488367203,
			48378511622303, 50577534878987, 52776558134423, 54975581392583, 57174604644503, 59373627900407,
			61572651156383, 63771674412287, 65970697666967, 68169720924167, 74766790688867, 79164837200927,
			83562883712027, 87960930223163, 92358976733483, 96757023247427, 101155069756823, 105553116266999,
			109951162779203, 114349209290003, 118747255800179, 123145302311783, 127543348823027, 131941395333479,
			136339441846019, 149533581378263, 158329674402959, 167125767424739, 175921860444599, 184717953466703,
			193514046490343, 202310139514283, 211106232536699, 219902325558107, 228698418578879, 237494511600287,
			246290604623279, 255086697645023, 263882790666959, 272678883689987, 299067162755363, 316659348799919,
			334251534845303, 351843720890723, 369435906934019, 387028092977819, 404620279022447, 422212465067447,
			439804651111103, 457396837157483, 474989023199423, 492581209246163, 510173395291199, 527765581341227,
			545357767379483, 598134325510343, 633318697599023, 668503069688723, 703687441776707, 738871813866287,
			774056185954967, 809240558043419, 844424930134187, 879609302222207, 914793674313899, 949978046398607,
			985162418489267, 1020346790579903, 1055531162666507, 1090715534754863 }; //list of 450 sorted prime numbers;
	vector<UINT64> primeNumbers(array, array + sizeof array / sizeof array[0]);
	for (size_t i = 0; i < primeNumbers.size(); i++)
		if (primeNumbers.at(i) > number)
			return primeNumbers.at(i); // Return the smallest prime in the list larger than number.
	return number + 1;
}

UINT64 HashTableLongKey::hashFunction(const string & _subString) {
	UINT64 sum1 = 1, sum2 = 1, length = _subString.length(); // initial value of 1 to avoid zero when doing the multiplication.
	for (UINT64 i = 0; i < length; i++) { // We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
		if (i < 32)
			// sum1 is for the first 32 bp. bit shifted to left 2 bits.
			// Change the character to integer. A=Ox41=01000001
			//                                  C=0x43=01000011
			//                                  G=0x47=01000111
			//                                  T=0x54=01010100
			// Then, shift to right way 1 bit.
			// Then, bit and operation with 00000011
			// Then, it just have A=00, C=01, G=11,T=10
			sum1 = (sum1 << 2) | ((((UINT64) (_subString[i])) >> 1) & 3ull);
		else
			sum2 = (sum2 << 2) | ((((UINT64) (_subString[i])) >> 1) & 3ull);
	}
	// multiply two mod results, and mod again.
	return ((sum1 % hashTableSize) * (sum2 % hashTableSize)) % hashTableSize; // Modulus operation to get the index in the hash table.
}

UINT64 HashTableLongKey::hashFunction(const string & _sequence, INT64 _start, INT64 _length) {
	UINT64 sum1 = 1, sum2 = 1, length = _length; // initial value of 1 to avoid zero when doing the multiplication.
	for (UINT64 i = 0; i < length; i++) { // We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
		if (i < 32)
			// sum1 is for the first 32 bp. bit shifted to left 2 bits.
			// Change the character to integer. A=Ox41=01000001
			//                                  C=0x43=01000011
			//                                  G=0x47=01000111
			//                                  T=0x54=01010100
			// Then, shift to right way 1 bit.
			// Then, bit and operation with 00000011
			// Then, it just have A=00, C=01, G=11,T=10
			sum1 = (sum1 << 2) | ((((UINT64) (_sequence[i + _start]) >> 1)) & 3ull);
		else
			sum2 = (sum2 << 2) | ((((UINT64) (_sequence[i + _start]) >> 1)) & 3ull);
	}
	// multiply two mod results, and mod again.
	return ((sum1 % hashTableSize) * (sum2 % hashTableSize)) % hashTableSize; // Modulus operation to get the index in the hash table.
}
/*

 UINT64 HashTableLongKey::hashFunction(const string & _subString) {
 UINT64 sum1 = 1, sum2 = 1, length = _subString.length(); // initial value of 1 to avoid zero when doing the multiplication.
 for (UINT64 i = 0; i < length; i++) { // We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
 if (i < 32)
 // sum1 is for the first 32 bp. bit shifted to left 2 bits.
 // Change the character to integer. A=Ox41=01000001
 //                                  C=0x43=01000011
 //                                  G=0x47=01000111
 //                                  T=0x54=01010100
 // Then, shift to right way 1 bit.
 // Then, bit and operation with 00000011
 // Then, it just have A=00, C=01, G=11,T=10
 sum1 = (sum1 << 2) | ((((UINT64) (_subString[i])) >> 1) & 3ull);
 else
 sum2 = (sum2 << 2) | ((((UINT64) (_subString[i])) >> 1) & 3ull);
 }
 // multiply two mod results, and mod again.
 return ((sum1 % hashTableSize) * (sum2 % hashTableSize)) % hashTableSize; // Modulus operation to get the index in the hash table.
 }

 UINT64 HashTableLongKey::hashFunction(const string & _sequence, INT64 _start, INT64 _length) {
 UINT64 sum1 = 1, sum2 = 1, length = _length; // initial value of 1 to avoid zero when doing the multiplication.
 for (UINT64 i = 0; i < length; i++) { // We take the bit representation of the string. A = 00, C = 01, G = 11 and T = 10
 if (i < 32)
 // sum1 is for the first 32 bp. bit shifted to left 2 bits.
 // Change the character to integer. A=Ox41=01000001
 //                                  C=0x43=01000011
 //                                  G=0x47=01000111
 //                                  T=0x54=01010100
 // Then, shift to right way 1 bit.
 // Then, bit and operation with 00000011
 // Then, it just have A=00, C=01, G=11,T=10
 sum1 = (sum1 << 2) | ((((UINT64) (_sequence[i + _start]) >> 1)) & 3ull);
 else
 sum2 = (sum2 << 2) | ((((UINT64) (_sequence[i + _start]) >> 1)) & 3ull);
 }
 // multiply two mod results, and mod again.
 return ((sum1 % hashTableSize) * (sum2 % hashTableSize)) % hashTableSize; // Modulus operation to get the index in the hash table.
 }
 */

bool HashTableLongKey::insertQueryDataset(QueryDataset* _qDataset) {

	CLOCKSTART
	MEMORYSTART
	int totalKey = 0;
	int singleKey = 0;
	int tempSingleKey = 0;
	if (this->vviHashTable == NULL) {
		cout << "Hash Table hasn't been created yet." << endl;
		return false;
	} else {
		UINT64 datasetsize = this->dataSet->getNumberOfReads();
		UINT64 currentID = 0;
		while (currentID < datasetsize) {
			if ((currentID + 1L) % 1000000L == 0)
				cout << "Number of reads inserted in the hash table: " << (currentID + 1L) << "\r";
			QueryRead * read = this->dataSet->getReadFromID(currentID);
			if (read == NULL) {
				cout << endl << "Empty read: " << currentID << " ";
				currentID++;
				continue;
			}
			if (!read->bIsGoodRead) {
				currentID++;
				continue;
			}
			string * forwardRead = read->getSequence();
			INT64 i = 0;
			tempSingleKey = 0;
			for (i = 0; i <= read->getRawReadLength() - this->hashKeyLength; i += this->hashKeyLength) {
				string subString = forwardRead->substr(i, this->hashKeyLength);
				this->insertQueryRead(read, subString, i);
				tempSingleKey++;
			}
			if (i < read->getRawReadLength()) {
				i = read->getRawReadLength() - this->hashKeyLength;
				string subString = forwardRead->substr(i, this->hashKeyLength);
				this->insertQueryRead(read, subString, i);
				tempSingleKey++;
			}
			currentID++;
			if (tempSingleKey > singleKey)
				singleKey = tempSingleKey;
			totalKey += tempSingleKey;
		}
		cout << endl;
		cout << "Hash Table" << " maximum collision number is: " << this->numberOfHashCollision << endl;
		cout << "Hash Table" << " maximum single read collision number is: " << this->maxSingleHashCollision << endl;
		cout << "Total number of Keys: " << totalKey << endl;
		cout << "Maximum number of Keys for a read " << singleKey << endl;

		MEMORYSTOP
		CLOCKSTOP

		return true;
	}
}

bool HashTableLongKey::insertQueryRead(QueryRead *read, string subString, INT64 postion) {
	UINT64 currentCollision = 0;

	UINT64 index = this->hashFunction(subString);
	while (this->vviHashTable->at(index) != NULL && (!this->vviHashTable->at(index)->empty())) {
		vector<INT64>* readList = this->vviHashTable->at(index);
		INT64 data = readList->at(0);
		INT64 keyreadID = data & this->iMask;
		INT64 keymode = (data >> this->iBitOfReadId);
		QueryRead * read = this->dataSet->getReadFromID(keyreadID);
		string keyStr = read->getSequence()->substr(keymode, this->hashKeyLength);
		if (strcmp(keyStr.c_str(), subString.c_str()) == 0)
			break;
		numberOfHashCollision++;
		currentCollision++;
		index = (index == this->hashTableSize - 1) ? 0 : index + 1; // Increment the index
		if (currentCollision >= this->hashTableSize) {
			cout << "The size of hash table is too small." << endl;
			exit(1);
		}
	}
	if (this->vviHashTable->at(index) == NULL) {
		vector<INT64> * newList = new vector<INT64>;
		newList->resize(newList->size());
		this->vviHashTable->at(index) = newList;
	}
	INT64 positionShifted = (postion << this->iBitOfReadId);
	INT64 iReadInnerId = read->getReadInnerId();
	INT64 key = positionShifted | iReadInnerId;
	this->vviHashTable->at(index)->push_back(key); // Add the string in the list.

	if (currentCollision > this->maxSingleHashCollision)
		this->maxSingleHashCollision = currentCollision;
	return true;
}

bool HashTableLongKey::isEmptyAt(INT64 _hashTableIndex) {
	if (this->vviHashTable->at(_hashTableIndex) == NULL) {
		return true;
	}
	return this->vviHashTable->at(_hashTableIndex)->empty();
}

void HashTableLongKey::setHashTableSizeAndInitialize(INT64 size) {
	cout << "Hash Table size set to: " << size << endl;

	this->vviHashTable = new vector<vector<INT64>*>(size, NULL);

	/*this->vviHashTable->reserve(size);

	 for (UINT64 i = 0; i < this->vviHashTable->capacity(); i++) { // Initialize the hash table.
	 vector<UINT64> * newList = new vector<UINT64>;
	 newList->resize(newList->size());
	 this->vviHashTable->push_back(newList);
	 }*/
}

