#ifndef BISON_DYN_PARSER_H
# define BISON_DYN_PARSER_H

#ifndef YYSTYPE
typedef union {
  double real;
  char* string;
  dyn* dyn_ptr;
  vector<double>* realvec_ptr;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	PARTICLE	257
# define	LOG	258
# define	DYNAMICS	259
# define	HYDRO	260
# define	STAR	261
# define	KEYWORD	262
# define	LOG_STORY	263
# define	STRING	264
# define	NUMBER	265


extern YYSTYPE yylval;

#endif /* not BISON_DYN_PARSER_H */
