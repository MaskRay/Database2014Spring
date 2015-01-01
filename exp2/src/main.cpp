#include "SimJoiner.h"
#include <cassert>
#include <unistd.h>
#include <err.h>
#include <cstdio>
#include <cstring>

int main(int argc, char **argv)
{
  if (argc != 6) {
    fprintf(stderr, "Usage: %s [ed|j] q threshold file1 file2\n", argv[0]);
    return 1;
  }
  if (access(argv[5], R_OK))
    errx(2, "%s not readable", argv[5]);

  size_t len = 0;
  char *line = NULL;
  VS doc1, doc2;

  FILE *fp = fopen(argv[4], "r");
  if (! fp)
    errx(2, "%s not readable", argv[4]);
  while (getline(&line, &len, fp) != -1) {
    line[strlen(line)-1] = '\0';
    doc1.pb(line);
  }

  fp = fopen(argv[5], "r");
  if (! fp)
    errx(2, "%s not readable", argv[5]);
  while (getline(&line, &len, fp) != -1) {
    line[strlen(line)-1] = '\0';
    doc2.pb(line);
  }

  bool debug = getenv("D") != NULL;

  SimJoiner joiner;
  string op = argv[1];
  if (op == "j") {
    vector<JaccardJoinResult> res;
    u q = atoi(argv[2]);
    ft th = atof(argv[3]);
    joiner.joinJaccard(argv[4], argv[5], q, th, res);

    if (debug)
      for (auto i: res)
        printf("+ %lf\n  %s\n  %s\n", i.s, doc1[i.id1].c_str(), doc2[i.id2].c_str());
    else
      for (auto i: res)
        printf("+ %u %u %lf\n", i.id1, i.id2, i.s);
  } else if (op == "ed") {
    vector<EDJoinResult> res;
    u q = atoi(argv[2]);
    u th = atoi(argv[3]);
    joiner.joinED(argv[4], argv[5], q, th, res);

    if (debug)
      for (auto i: res)
        printf("+ %u\n  %s\n  %s\n", i.s, doc1[i.id1].c_str(), doc2[i.id2].c_str());
    else
      for (auto i: res)
        printf("+ %u %u %u\n", i.id1, i.id2, i.s);
  } else
    errx(3, "unknown operation");

  free(line);

  return 0;
}
