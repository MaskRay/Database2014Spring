#include "SimSearcher.h"
#include <cassert>
#include <cstdio>
#include <cstring>

int main(int argc, char *argv[])
{
  char buf[2050];
  if (argc != 4) {
    fprintf(stderr, "Usage: %s q data query\n", argv[0]);
    return 1;
  }
  SimSearcher searcher;
  int q = atoi(argv[1]);
  vector<string> strs;
	vector<pair<unsigned, double> > res;

  searcher.createIndex(argv[2], q);
  FILE *fp = fopen(argv[2], "r");
  char *line = NULL;
  size_t len = 0;
  double threshold;
  while (getline(&line, &len, fp) != -1) {
    line[strlen(line)-1] = '\0';
    strs.push_back(line);
  }
  fclose(fp);

  fp = fopen(argv[3], "r");
  while (getline(&line, &len, fp) != -1) {
    if (sscanf(line, "%[^\t]\t%lf", buf, &threshold) != 2)
      return 2;
    searcher.searchJaccard(buf, threshold, res);
    printf("+ %lf %s\n", threshold, buf);
    for (auto s: res)
      printf("  %lf %s\n", s.second, strs[s.first].c_str());
  }
  fclose(fp);
}
