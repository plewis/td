#pragma once

//#define DEBUGGING
//#define CLUSTER_DISTANCE
//#define KUHNER_FELSENSTEIN_DISTANCE
//#define DEBUGGING_LAPJV
//#define TESTKDE
//#define OP_SAVE_DOT_FILE

// This was introduced to fix the unrooted bug in which BHV distances
// calculated for unrooted trees failed to account for the root tip edge.
// The solution enabled by ALWAYS_ROOTED is to simply always treat trees
// as rooted. If ALWAYS_ROOTED is undefined, the unrooted bug will return!
#define ALWAYS_ROOTED
