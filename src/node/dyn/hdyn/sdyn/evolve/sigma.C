
// sigma.C: outline only...

#include "sdyn.h"

typedef struct _element {char* id;
			 int counter;
			 struct _element *next;	// singly-linked list
		     } list_element;

list_element* locate_string(list_element* e, char* string)
{
    // Find the list element indexed by the specified string.

    while (e) {
	if (!strcmp(string, e->id)) return e;
	e = e->next;
    }
    return NULL;
}

list_element* add_element(list_element* e, char* string)
{
    // Add a new element to the end of the list.

    if (e) {
	while (e->next) e = e->next;
	list_element* add = new list_element;
	e->next = add;
	e = add;
    } else
	e = new list_element;
    
    e->id = new char[strlen(string)];
    strcpy(e->id, string);

    e->counter = 1;
    e->next = NULL;

    return e;
}

void print_results(list_element* e)
{
    while (e) {
	cerr << e->id << "  " << e->counter << endl;
	e = e->next;
    }
}

main()
{
    list_element* start = NULL;

    for (;;) { 	// loop over experiments

	// Perform the experiment:

	sdyn* b = next_system(...);
	scatter(b,...);

	// Determine the normal form of the outcome:

	char* normal_form = get_normal_form(b);

	// Add the results to the linked list:

	list_element* e = locate_string(start, normal_form);

	if (!e) {

	    // Add a new element, or initialize the list:

	    list_element* temp = add_element(start, normal_form);
	    if (start == NULL) start = temp;

	} else

	    e->counter++;	// Just increment the counter
    }

    print_results(start);
}
