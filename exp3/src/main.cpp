#include "AEE.h"

using namespace std;

int main(int argc, char **argv)
{
	AEE aee;

	vector<EDExtractResult> resultED;
	vector<JaccardExtractResult> resultJaccard;

  bool debug = getenv("D") || getenv("V");

  if (argc != 4) {
    fprintf(stderr, "%s threshold e d\n", argv[0]);
    return 1;
  }

	unsigned q = 3, edThreshold = atoi(argv[1]);
	double jaccardThreshold = 0.85;

	aee.createIndex(argv[2], q);
	//aee.aeeJaccard(argv[2], jaccardThreshold, resultJaccard);

  char line[15000];
  auto tup = read_file(argv[3]);
  char *mapped = get<1>(tup);
  int file_size = int(get<2>(tup)), did = 0;
  for (int i = 0, j = 0; j < file_size; i = j, did++) {
    for (; j < file_size && mapped[j] != '\n'; j++);
    memcpy(line, mapped+i, j-i);
    line[j-i] = 0;
    if (debug)
      printf("+ %s\n", argv[2]);
    if (argv[3] && argv[3][0] == 'j') {
      aee.aeeJaccard(line, jaccardThreshold, resultJaccard);
      for (auto &r: resultJaccard)
        printf("%d %d %d %d %lf\n", did, r.id, r.pos, r.len, r.sim);
    } else {
      aee.aeeED(line, edThreshold, resultED);
      for (auto &r: resultED)
        printf("%d %d %d %d %d\n", did, r.id, r.pos, r.len, r.sim);
    }
    if (j < file_size && mapped[j] == '\n')
      j++;
  }

	return 0;
}

