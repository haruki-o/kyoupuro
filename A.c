#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Domain {
  char name[256];
  struct Domain *next;
};

struct Domain *list_add(struct Domain *list, char *t) {
  /* 新しい要素を作成 */
  struct Domain *node = (struct Domain *)malloc(sizeof(struct Domain));
  strcpy(node->name, t);
  node->next = NULL;

  if (list == NULL) {
    node->next = NULL;
    return node;
  } else {
    /* リストの末尾の要素を探す */
    struct Domain *p = list;
    while (p->next != NULL) {
      p = p->next;
    }
    /* ここで変数 p は末尾の要素をさす */
    p->next = node;
    return list;
  }
}

void list_show(struct Domain *list) {
  if (list != NULL) {
    /* リストの末尾の要素を探す */
    struct Domain *p = list;
    printf("{");
    while (1) {
      printf("%s, ", p->name);
      if (p->next == NULL)
        break;
      p = p->next;
    }
    printf("}\n");
  }
}

struct Domain *list_add_init(struct Domain *list, int argc, char *argv[]) {
  if (argc == 1) {

  } else {
    for (int i = 1; i < argc; i++) {
      FILE *fp;
      char str[256];
      fp = fopen(argv[i], "r");

      if (fp == NULL) {
        printf("ファイルオープンエラー");
        exit(0);
      }

      while (1) {
        // while ((fgets(str, sizeof(str), fp)) != NULL) {
        if ((fgets(str, sizeof(str), fp)) == NULL) {
          *str = '\0';
          break;
        }
        if (*str && str[strlen(str) - 1] == '\n') {
          str[strlen(str) - 1] = 0;
        }

        list = list_add(list, str);
        list_show(list);
      }
      fclose(fp);
    }
  }
  return list;
}

int list_size(struct Domain *list) {
  int ret = 0;
  if (list != NULL) {
    /* リストの末尾の要素を探す */
    struct Domain *p = list;
    while (1) {
      ret++;
      if (p->next == NULL)
        break;
      p = p->next;
    }
  }
  return ret;
}

int dot_count(char *s) {
  int ret = 0;
  for (int i = 0; i < strlen(s); i++) {
    if (s[i] == '.')
      ret++;
  }
  return ret;
}

int list_strcmp(char s1[], char s2[]) {
  int d1 = dot_count(s1);
  int d2 = dot_count(s2);

  char *domain_vec1[d1], *domain_vec2[d2];
  printf("%s, %s\n", s1, s2);
  int at_idx = d1;
  int l = 0;
  for (int i = 0; i < strlen(s1); i++) {
    if (s1[i] == '.') {
      substring(s1, l, i - l, domain_vec1[at_idx]);
      at_idx--;
      l = l + 1;
    }
  }

  // *at_s = '\0';
  // at_idx = d2 - 1;
  // for (int i = 0; i < strlen(s2); i++) {
  //   if (s2[i] == '.') {
  //     domain_vec2[at_idx] = at_s;
  //     at_idx--;
  //     at_s = '\0';
  //   } else {
  //     at_s += s2[i];
  //   }
  // }

  printf("%s\n", s1);
  for (int i = 0; i < d1 + 1; i++) {
    printf("{%s}", domain_vec1[i]);
  }
  // printf("%s\n", s2);
  // for (int i = 0; i < d2 + 1; i++) {
  //   printf("{%s}", domain_vec2[i]);
  // }

  return 0;
}

void list_sort(struct Domain *list) {
  int sum = list_size(list);
  printf("%d\n", sum);
  struct Domain *ip = list;
  char *best_str = ip->name;
  for (int i = 0; i < sum; i++) {
    // [i,sum)での辞書順先頭
    struct Domain *jp = ip;
    jp = jp->next;
    for (int j = i + 1; j < sum; j++) {
      if (list_strcmp(best_str, jp->name)) {
      }
      jp = jp->next;
    }
    printf("%d : %s\n", i, best_str);
    ip = ip->next;
  }
}

int main(int argc, char *argv[]) {
  struct Domain *list = NULL;
  list = list_add_init(list, argc, argv);
  list_sort(list);
  return 0;
}
