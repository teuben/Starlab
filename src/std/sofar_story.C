//
//  story.C
//

#include "story.h"

story::~story()
    {
    story * si = first_daughter_node;
    story * sn;
    while (si != NULL)
	{
	sn = si->next_story_node;
	delete si;
	si = sn;
	}
    delete [] text;
    }

int  is_chapter(story * s)
    {
    return s->get_chapter_flag();
    }

int  is_chapter_begin_line(char * line)
    {
    return (*line == chapter_begin_char);
    }

int  is_chapter_end_line(char * line)
    {
    return (*line == chapter_end_char);
    }

int get_story_line(istream & str, char * line)
    {
    str.get(line,MAX_STORY_LINE_LENGTH,'\n');

    if(str.eof())
        return 0;

    char c;
    if(str.get(c) && c!='\n') 
	{
	cerr << "get_story_line : input line too long :'"<<line<<"'\n";
	exit(1);
	}

    return 1;
    }

void add_daughter_story(story * s, story * d)
    {
    story *lc = s->get_last_daughter_node();
    if (lc)
	{
	lc->set_next_story_node(d);
	s->set_last_daughter_node(d);
	}
    else
	{
	s->set_first_daughter_node(d);
	s->set_last_daughter_node(d);
	}
    }

void add_chapter(story * s, story * chap)
    {
    if (chap == NULL)
	return;

    add_daughter_story(s, chap);
    }

void rm_daughter_story(story * s, story * d)
    {
    story * fd = s->get_first_daughter_node();
    story * ad;                                    // a daughter
    story * nd;                                    // next daughter

    if (fd == d)
	s->set_first_daughter_node(d->get_next_story_node());
    else
	{
	ad = fd;
	nd = ad->get_next_story_node();
	while (nd != d)
	    {
	    ad = nd;
	    nd = ad->get_next_story_node();
	    }
	nd = d->get_next_story_node();
	ad->set_next_story_node(nd);
	if (nd == NULL)
	    s->set_last_daughter_node(ad);
	}
    delete d;    
    }

story* mk_story_line()    {story* s = new story(0); return s;}

story* mk_story_line(char * line)
    {
    story* s = new story(0);
    s->set_text(line);
    return s;
    }

story* mk_story_chapter()       {story* s = new story(1); return s;}

story* mk_story_chapter(char * title)
    {
    story* s = new story(1);
    s->set_text(title);
    return s;
    }

void add_story_line(story * s, char * line)
    {
    if (line == NULL)
	return;

    story * new_s = mk_story_line(line);
    add_daughter_story(s, new_s);
    }

story* get_chapter(istream& str, char* line)
    {
    if (!is_chapter_begin_line(line))
        {
	cerr << "get_chapter: first line not chapter_begin_line";
	exit(1);
	}

    story * chap = mk_story_chapter(++line);

    char new_line[MAX_STORY_LINE_LENGTH];

    while (get_story_line(str, new_line) && !is_chapter_end_line(new_line))
	{
	if (is_chapter_begin_line(new_line))
	    add_chapter(chap, get_chapter(str, new_line));
	else
	    add_story_line(chap, new_line);
	}

    if (new_line == NULL)
	{
	cerr << "get_chapter: new_line == NULL before end of chapter\n";
	exit(1);
	}

    if (!streq(new_line+1, chap->get_text()))
	{
	cerr << "get_chapter: closing title ``" << new_line+1
	     << "'' differs from opening title ``" << chap->get_text()
	     << "''\n";
	exit(1);
	}

    return chap;
    }

story* get_story(istream & str, char *line)  {return get_chapter(str, line);}

story* get_story(istream& str)
    { 
    char line[MAX_STORY_LINE_LENGTH];

    if (!get_story_line(str, line))
        return(NULL);

    if (!is_chapter_begin_line(line))
        {
	cerr << "get_story: first line not chapter_begin_line\n";
	exit(1);
	}

    return get_chapter(str, line);
    }

void put_headline(ostream& str, story& s)
    {
    str << chapter_begin_char << s.get_text() << "\n";
    }

void put_tailline(ostream& str, story& s)
    {
    str << chapter_end_char << s.get_text() << "\n";
    }

void put_line_text(ostream& str, story& s)
    {
    str << s.get_text() << "\n";
    }

void put_chapter(ostream& str, story& s)
    {
    if (!is_chapter(&s))
	{
	cerr << "put_chapter: not a story\n";
	exit(1);
	}
    put_headline(str, s);

    for (story * d = s.get_first_daughter_node(); d != NULL;
						  d = d->get_next_story_node())
	{
	if (is_chapter(d))
	    put_chapter(str, *d);
	else
	    put_line_text(str, *d);
	}

    put_tailline(str, s);
    }

void put_story_contents(ostream& str, story& s)
    {
    if (!is_chapter(&s))
	{
	cerr << "put_story_contents: not a story\n";
	exit(1);
	}

    for (story * d = s.get_first_daughter_node(); d != NULL;
						  d = d->get_next_story_node())
	{
	if (is_chapter(d))
	    put_chapter(str, *d);
	else
	    put_line_text(str, *d);
	}
    }

void put_story(ostream& str, story& s)
    {
    put_chapter(str, s);
    }

ostream & operator << (ostream & s, story * sptr)
    {put_story(s, *sptr); return s;}

istream & operator >> (istream & s, story * & sptr)
    {sptr = get_story(s); return s;}

#ifdef TOOLBOX

main()
    {
    story * s;

    while (cin >> s)
	{
#if 0
	putiq(s, "test message", 41);
	putiq(s, "test message", 42);
	putiq(s, "yet another test message", 137);
	putrq(s, "pi", 3.14);
	putsq(s, "star cluster", "47 Tuc");
	vec * tmp = new vec(1,2,3);
	putvq(s, "1-2-3-test vector", *tmp);
#endif
//	cout << "test message = " << getiq(s, "test message") << endl;
//	cout << "yet another test message = "
//	     << getiq(s, "yet another test message") << endl;
//	cout << "pi = " << getrq(s, "pi") << endl;
//	cout << "star cluster = " << getsq(s, "star cluster") << endl;
//	cout << "1-2-3-test vector = " << getvq(s, "1-2-3-test vector")
//	     << endl;
        cout << s;
	}

//  cerr << "story.C : TOOLBOX : Normal exit\n";
    }

#endif
