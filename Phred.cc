#include "Phred.hh"


// Generated with this Python:
// s='{'
// for i in xrange(0,33) :
//     s+= '42,'
// for i in xrange(33,75) :
//     s+= str(i-33) + ','
// for i in xrange(75,127) :
//     s+= '42,'
// s+='42}'
const size_t Phred::decode_phred33_lut[128] = {
    42,42,42,42,42,42,42,42,42,42,
    42,42,42,42,42,42,42,42,42,42,
    42,42,42,42,42,42,42,42,42,42,
    42,42,42, 0, 1, 2, 3, 4, 5, 6,
    7 ,8 ,9 ,10,11,12,13,14,15,16,
    17,18,19,20,21,22,23,24,25,26,
    27,28,29,30,31,32,33,34,35,36,
    37,38,39,40,41,42,42,42,42,42,
    42,42,42,42,42,42,42,42,42,42,
    42,42,42,42,42,42,42,42,42,42,
    42,42,42,42,42,42,42,42,42,42,
    42,42,42,42,42,42,42,42,42,42,
    42,42,42,42,42,42,42,42
};

// Generated with this Python:
// ','.join(('"'+chr(i)+'"' for i in xrange(33,33+41)))
const char* Phred::encode_phred33_lut[NUM_QUALS+1] = {
    "'!'","'\"'","'#'","'$'","'%'","'&'","'\\''","'('","')'","'*'",
    "'+'","','","'-'","'.'","'/'","'0'","'1'","'2'","'3'","'4'",
    "'5'","'6'","'7'","'8'","'9'","':'","';'","'<'","'='","'>'",
    "'?'","'@'","'A'","'B'","'C'","'D'","'E'","'F'","'G'","'H'",
    "'I'", "'J'", "NA"
};
