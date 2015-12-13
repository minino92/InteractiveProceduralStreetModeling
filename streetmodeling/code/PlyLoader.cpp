#include "stdafx.h"
#include ".\plyloader.h"

char *type_names[] = {  /* names of scalar types */
"invalid",
"int8", "int16", "int32", "uint8", "uint16", "uint32", "float32", "float64",
};

char *old_type_names[] = {  /* old names of types for backward compatability */
"invalid",
"char", "short", "int", "uchar", "ushort", "uint", "float", "double",
};

int ply_type_size[] = {
  0, 1, 2, 4, 1, 2, 4, 4, 8
};

static char *my_alloc(int size, int lnum, char *fname)
{
  char *ptr;

  ptr = (char *) malloc (size);

  if (ptr == 0) {
    fprintf(stderr, "Memory allocation bombed on line %d in %s\n", lnum, fname);
  }

  return (ptr);
}


/* memory allocation */

#define myalloc(mem_size) my_alloc((mem_size), __LINE__, __FILE__)



CPlyLoader::CPlyLoader(void)
{
}

CPlyLoader::~CPlyLoader(void)
{
}

char *CPlyLoader::setup_element_read_ply (PlyFile *ply, int index, int *elem_count)
{
  PlyElement *elem;

  if (index < 0 || index > ply->num_elem_types) {
    fprintf (stderr, "Warning:  No element with index %d\n", index);
    return (0);
  }

  elem = ply->elems[index];

  /* set this to be the current element */
  ply->which_elem = elem;

  /* return the number of such elements in the file and the element's name */
  *elem_count = elem->num;
  return (elem->name);
}

int CPlyLoader::equal_strings(char *s1, char *s2)
{

  while (*s1 && *s2)
    if (*s1++ != *s2++)
      return (0);

  if (*s1 != *s2)
    return (0);
  else
    return (1);
}


PlyProperty *CPlyLoader::find_property(PlyElement *elem, char *prop_name, int *index)
{
  int i;

  for (i = 0; i < elem->nprops; i++)
    if (equal_strings (prop_name, elem->props[i]->name)) {
      *index = i;
      return (elem->props[i]);
    }

  *index = -1;
  return (NULL);
}


void CPlyLoader::setup_property_ply(
  PlyFile *plyfile,
  PlyProperty *prop
  )
{
  PlyElement *elem;
  PlyProperty *prop_ptr;
  int index;

  elem = plyfile->which_elem;

  /* deposit the property information into the element's description */

  prop_ptr = find_property (elem, prop->name, &index);
  if (prop_ptr == NULL) {
    fprintf (stderr, "Warning:  Can't find property '%s' in element '%s'\n",
             prop->name, elem->name);
    return;
  }
  prop_ptr->internal_type  = prop->internal_type;
  prop_ptr->offset         = prop->offset;
  prop_ptr->count_internal = prop->count_internal;
  prop_ptr->count_offset   = prop->count_offset;

  /* specify that the user wants this property */
  elem->store_prop[index] = STORE_PROP;

}

void CPlyLoader::setup_other_props(PlyFile *plyfile, PlyElement *elem)
{
  int i;
  PlyProperty *prop;
  int size = 0;
  int type_size;

  /* Examine each property in decreasing order of size. */
  /* We do this so that all data types will be aligned by */
  /* word, half-word, or whatever within the structure. */

  for (type_size = 8; type_size > 0; type_size /= 2) {

    /* add up the space taken by each property, and save this information */
    /* away in the property descriptor */

    for (i = 0; i < elem->nprops; i++) {

      /* don't bother with properties we've been asked to store explicitly */
      if (elem->store_prop[i])
        continue;

      prop = elem->props[i];

      /* internal types will be same as external */
      prop->internal_type = prop->external_type;
      prop->count_internal = prop->count_external;

      /* list case */
      if (prop->is_list == PLY_LIST) {

        /* pointer to list */
        if (type_size == sizeof (void *)) {
          prop->offset = size;
          size += sizeof (void *);    /* always use size of a pointer here */
        }

        /* count of number of list elements */
        if (type_size == ply_type_size[prop->count_external]) {
          prop->count_offset = size;
          size += ply_type_size[prop->count_external];
        }
      }
      /* string */
      else if (prop->is_list == PLY_STRING) {
        /* pointer to string */
        if (type_size == sizeof (char *)) {
          prop->offset = size;
          size += sizeof (char *);
        }
      }
      /* scalar */
      else if (type_size == ply_type_size[prop->external_type]) {
        prop->offset = size;
        size += ply_type_size[prop->external_type];
      }
    }

  }

  /* save the size for the other_props structure */
  elem->other_size = size;
}


PlyOtherProp *CPlyLoader::get_other_properties_ply(
  PlyFile *plyfile,
  int offset
)
{
  PlyOtherProp *other;

  other = get_other_properties (plyfile, plyfile->which_elem, offset);
  return (other);
}

 PlyOtherProp *CPlyLoader::get_other_properties(
  PlyFile *plyfile,
  PlyElement *elem,
  int offset
)
{
  int i;
  PlyOtherProp *other;
  PlyProperty *prop;
  int nprops;

  /* remember that this is the "current" element */
  plyfile->which_elem = elem;

  /* save the offset to where to store the other_props */
  elem->other_offset = offset;

  /* place the appropriate pointers, etc. in the element's property list */
  setup_other_props (plyfile, elem);

  /* create structure for describing other_props */
  other = (PlyOtherProp *) myalloc (sizeof (PlyOtherProp));
  other->name = strdup (elem->name);
#if 0
  if (elem->other_offset == NO_OTHER_PROPS) {
    other->size = 0;
    other->props = NULL;
    other->nprops = 0;
    return (other);
  }
#endif
  other->size = elem->other_size;
  other->props = (PlyProperty **) myalloc (sizeof(PlyProperty) * elem->nprops);
  
  /* save descriptions of each "other" property */
  nprops = 0;
  for (i = 0; i < elem->nprops; i++) {
    if (elem->store_prop[i])
      continue;
    prop = (PlyProperty *) myalloc (sizeof (PlyProperty));
    copy_property (prop, elem->props[i]);
    other->props[nprops] = prop;
    nprops++;
  }
  other->nprops = nprops;

  /* set other_offset pointer appropriately if there are NO other properties */
  if (other->nprops == 0) {
    elem->other_offset = NO_OTHER_PROPS;
  }
 
  /* return structure */
  return (other);
}


void CPlyLoader::get_element_ply (PlyFile *plyfile, void *elem_ptr)
{
  if (plyfile->file_type == PLY_ASCII)
    ascii_get_element (plyfile, (char *) elem_ptr);
  else
    binary_get_element (plyfile, (char *) elem_ptr);
}

void CPlyLoader::ascii_get_element(PlyFile *plyfile, char *elem_ptr)
{
  int i,j,k;
  PlyElement *elem;
  PlyProperty *prop;
  char **words;
  int nwords;
  int which_word;
  FILE *fp = plyfile->fp;
  char *elem_data,*item;
  char *item_ptr;
  int item_size;
  int int_val;
  unsigned int uint_val;
  double double_val;
  int list_count;
  int store_it;
  char **store_array;
  char *orig_line;
  char *other_data;
  int other_flag;

  /* the kind of element we're reading currently */
  elem = plyfile->which_elem;

  /* do we need to setup for other_props? */

  if (elem->other_offset != NO_OTHER_PROPS) {
    char **ptr;
    other_flag = 1;
    /* make room for other_props */
    other_data = (char *) myalloc (elem->other_size);
    /* store pointer in user's structure to the other_props */
    ptr = (char **) (elem_ptr + elem->other_offset);
    *ptr = other_data;
  }
  else
    other_flag = 0;

  /* read in the element */

  words = get_words (plyfile->fp, &nwords, &orig_line);
  if (words == NULL) {
    fprintf (stderr, "ply_get_element: unexpected end of file\n");
    exit (-1);
  }

  which_word = 0;

  for (j = 0; j < elem->nprops; j++) {

    prop = elem->props[j];
    store_it = (elem->store_prop[j] | other_flag);

    /* store either in the user's structure or in other_props */
    if (elem->store_prop[j])
      elem_data = elem_ptr;
    else
      elem_data = other_data;

    if (prop->is_list == PLY_LIST) {       /* a list */

      /* get and store the number of items in the list */
      get_ascii_item (words[which_word++], prop->count_external,
                      &int_val, &uint_val, &double_val);
      if (store_it) {
        item = elem_data + prop->count_offset;
        store_item(item, prop->count_internal, int_val, uint_val, double_val);
      }

      /* allocate space for an array of items and store a ptr to the array */
      list_count = int_val;
      item_size = ply_type_size[prop->internal_type];
      store_array = (char **) (elem_data + prop->offset);

      if (list_count == 0) {
        if (store_it)
          *store_array = NULL;
      }
      else {
        if (store_it) {
          item_ptr = (char *) myalloc (sizeof (char) * item_size * list_count);
          item = item_ptr;
          *store_array = item_ptr;
        }

        /* read items and store them into the array */
        for (k = 0; k < list_count; k++) {
          get_ascii_item (words[which_word++], prop->external_type,
                          &int_val, &uint_val, &double_val);
          if (store_it) {
            store_item (item, prop->internal_type,
                        int_val, uint_val, double_val);
            item += item_size;
          }
        }
      }

    }
    else if (prop->is_list == PLY_STRING) {   /* a string */
      if (store_it) {
	char *str;
	char **str_ptr;
	str = strdup (words[which_word++]);
        item = elem_data + prop->offset;
	str_ptr = (char **) item;
	*str_ptr = str;
      }
      else {
        which_word++;
      }
    }
    else {                     /* a scalar */
      get_ascii_item (words[which_word++], prop->external_type,
                      &int_val, &uint_val, &double_val);
      if (store_it) {
        item = elem_data + prop->offset;
        store_item (item, prop->internal_type, int_val, uint_val, double_val);
      }
    }

  }

  free (words);
}


char **CPlyLoader::get_words(FILE *fp, int *nwords, char **orig_line)
{
#define BIG_STRING 4096
  int i,j;
  static char str[BIG_STRING];
  static char str_copy[BIG_STRING];
  char **words;
  int max_words = 10;
  int num_words = 0;
  char *ptr,*ptr2;
  char *result;

  words = (char **) myalloc (sizeof (char *) * max_words);

  /* read in a line */
  result = fgets (str, BIG_STRING, fp);
  if (result == NULL) {
    *nwords = 0;
    *orig_line = NULL;
    return (NULL);
  }

  /* convert line-feed and tabs into spaces */
  /* (this guarentees that there will be a space before the */
  /*  null character at the end of the string) */

  str[BIG_STRING-2] = ' ';
  str[BIG_STRING-1] = '\0';

  for (ptr = str, ptr2 = str_copy; *ptr != '\0'; ptr++, ptr2++) {
    *ptr2 = *ptr;
    if (*ptr == '\t') {
      *ptr = ' ';
      *ptr2 = ' ';
    }
    else if (*ptr == '\n') {
      *ptr = ' ';
      *ptr2 = '\0';
      break;
    }
  }

  /* find the words in the line */

  ptr = str;
  while (*ptr != '\0') {

    /* jump over leading spaces */
    while (*ptr == ' ')
      ptr++;

    /* break if we reach the end */
    if (*ptr == '\0')
      break;

    /* allocate more room for words if necessary */
    if (num_words >= max_words) {
      max_words += 10;
      words = (char **) realloc (words, sizeof (char *) * max_words);
    }

    if (*ptr == '\"') {  /* a quote indidicates that we have a string */

      /* skip over leading quote */
      ptr++;

      /* save pointer to beginning of word */
      words[num_words++] = ptr;

      /* find trailing quote or end of line */
      while (*ptr != '\"' && *ptr != '\0')
        ptr++;

      /* replace quote with a null character to mark the end of the word */
      /* if we are not already at the end of the line */
      if (*ptr != '\0')
	*ptr++ = '\0';
    }
    else {               /* non-string */

      /* save pointer to beginning of word */
      words[num_words++] = ptr;

      /* jump over non-spaces */
      while (*ptr != ' ')
	ptr++;

      /* place a null character here to mark the end of the word */
      *ptr++ = '\0';
    }
  }

  /* return the list of words */
  *nwords = num_words;
  *orig_line = str_copy;
  return (words);
}


void CPlyLoader::get_ascii_item(
  char *word,
  int type,
  int *int_val,
  unsigned int *uint_val,
  double *double_val
)
{
  switch (type) {
    case Int8:
    case Uint8:
    case Int16:
    case Uint16:
    case Int32:
      *int_val = atoi (word);
      *uint_val = *int_val;
      *double_val = *int_val;
      break;

    case Uint32:
      *uint_val = strtoul (word, (char **) NULL, 10);
      *int_val = *uint_val;
      *double_val = *uint_val;
      break;

    case Float32:
    case Float64:
      *double_val = atof (word);
      *int_val = (int) *double_val;
      *uint_val = (unsigned int) *double_val;
      break;

    default:
      fprintf (stderr, "get_ascii_item: bad type = %d\n", type);
      exit (-1);
  }
}

void CPlyLoader::store_item (
  char *item,
  int type,
  int int_val,
  unsigned int uint_val,
  double double_val
)
{
  unsigned char *puchar;
  short int *pshort;
  unsigned short int *pushort;
  int *pint;
  unsigned int *puint;
  float *pfloat;
  double *pdouble;

  switch (type) {
    case Int8:
      *item = int_val;
      break;
    case Uint8:
      puchar = (unsigned char *) item;
      *puchar = uint_val;
      break;
    case Int16:
      pshort = (short *) item;
      *pshort = int_val;
      break;
    case Uint16:
      pushort = (unsigned short *) item;
      *pushort = uint_val;
      break;
    case Int32:
      pint = (int *) item;
      *pint = int_val;
      break;
    case Uint32:
      puint = (unsigned int *) item;
      *puint = uint_val;
      break;
    case Float32:
      pfloat = (float *) item;
      *pfloat = double_val;
      break;
    case Float64:
      pdouble = (double *) item;
      *pdouble = double_val;
      break;
    default:
      fprintf (stderr, "store_item: bad type = %d\n", type);
      exit (-1);
  }
}


void  CPlyLoader::binary_get_element(PlyFile *plyfile, char *elem_ptr)
{
  int i,j,k;
  PlyElement *elem;
  PlyProperty *prop;
  FILE *fp = plyfile->fp;
  char *elem_data;
  char *item;
  char *item_ptr;
  int item_size;
  int int_val;
  unsigned int uint_val;
  double double_val;
  int list_count;
  int store_it;
  char **store_array;
  char *other_data;
  int other_flag;

  /* the kind of element we're reading currently */
  elem = plyfile->which_elem;

  /* do we need to setup for other_props? */

  if (elem->other_offset != NO_OTHER_PROPS) {
    char **ptr;
    other_flag = 1;
    /* make room for other_props */
    other_data = (char *) myalloc (elem->other_size);
    /* store pointer in user's structure to the other_props */
    ptr = (char **) (elem_ptr + elem->other_offset);
    *ptr = other_data;
  }
  else
    other_flag = 0;

  /* read in a number of elements */

  for (j = 0; j < elem->nprops; j++) {

    prop = elem->props[j];
    store_it = (elem->store_prop[j] | other_flag);

    /* store either in the user's structure or in other_props */
    if (elem->store_prop[j])
      elem_data = elem_ptr;
    else
      elem_data = other_data;

    if (prop->is_list == PLY_LIST) {          /* list */

      /* get and store the number of items in the list */
      get_binary_item (fp, prop->count_external,
                      &int_val, &uint_val, &double_val);
      if (store_it) {
        item = elem_data + prop->count_offset;
        store_item(item, prop->count_internal, int_val, uint_val, double_val);
      }

      /* allocate space for an array of items and store a ptr to the array */
      list_count = int_val;
      item_size = ply_type_size[prop->internal_type];
      store_array = (char **) (elem_data + prop->offset);
      if (list_count == 0) {
        if (store_it)
          *store_array = NULL;
      }
      else {
        if (store_it) {
          item_ptr = (char *) myalloc (sizeof (char) * item_size * list_count);
          item = item_ptr;
          *store_array = item_ptr;
        }

        /* read items and store them into the array */
        for (k = 0; k < list_count; k++) {
          get_binary_item (fp, prop->external_type,
                          &int_val, &uint_val, &double_val);
          if (store_it) {
            store_item (item, prop->internal_type,
                        int_val, uint_val, double_val);
            item += item_size;
          }
        }
      }

    }
    else if (prop->is_list == PLY_STRING) {     /* string */
      int len;
      char *str;
      fread (&len, sizeof(int), 1, fp);
      str = (char *) myalloc (len);
      fread (str, len, 1, fp);
      if (store_it) {
	char **str_ptr;
        item = elem_data + prop->offset;
	str_ptr = (char **) item;
	*str_ptr = str;
      }
    }
    else {                                      /* scalar */
      get_binary_item (fp, prop->external_type,
                      &int_val, &uint_val, &double_val);
      if (store_it) {
        item = elem_data + prop->offset;
        store_item (item, prop->internal_type, int_val, uint_val, double_val);
      }
    }

  }
}


void CPlyLoader::get_binary_item(
  FILE *fp,
  int type,
  int *int_val,
  unsigned int *uint_val,
  double *double_val
)
{
  char c[8];
  void *ptr;

  ptr = (void *) c;

  switch (type) {
    case Int8:
      fread (ptr, 1, 1, fp);
      *int_val = *((char *) ptr);
      *uint_val = *int_val;
      *double_val = *int_val;
      break;
    case Uint8:
      fread (ptr, 1, 1, fp);
      *uint_val = *((unsigned char *) ptr);
      *int_val = *uint_val;
      *double_val = *uint_val;
      break;
    case Int16:
      fread (ptr, 2, 1, fp);
      *int_val = *((short int *) ptr);
      *uint_val = *int_val;
      *double_val = *int_val;
      break;
    case Uint16:
      fread (ptr, 2, 1, fp);
      *uint_val = *((unsigned short int *) ptr);
      *int_val = *uint_val;
      *double_val = *uint_val;
      break;
    case Int32:
      fread (ptr, 4, 1, fp);
      *int_val = *((int *) ptr);
      *uint_val = *int_val;
      *double_val = *int_val;
      break;
    case Uint32:
      fread (ptr, 4, 1, fp);
      *uint_val = *((unsigned int *) ptr);
      *int_val = *uint_val;
      *double_val = *uint_val;
      break;
    case Float32:
      fread (ptr, 4, 1, fp);
      *double_val = *((float *) ptr);
      *int_val = *double_val;
      *uint_val = *double_val;
      break;
    case Float64:
      fread (ptr, 8, 1, fp);
      *double_val = *((double *) ptr);
      *int_val = *double_val;
      *uint_val = *double_val;
      break;
    default:
      fprintf (stderr, "get_binary_item: bad type = %d\n", type);
      exit (-1);
  }
}

PlyOtherElems *CPlyLoader::get_other_element_ply (PlyFile *plyfile)
{
  int i;
  PlyElement *elem;
  char *elem_name;
  int elem_count;
  PlyOtherElems *other_elems;
  OtherElem *other;

  elem = plyfile->which_elem;
  elem_name = elem->name;
  elem_count = elem->num;

  /* create room for the new "other" element, initializing the */
  /* other data structure if necessary */

  if (plyfile->other_elems == NULL) {
    plyfile->other_elems = (PlyOtherElems *) myalloc (sizeof (PlyOtherElems));
    other_elems = plyfile->other_elems;
    other_elems->other_list = (OtherElem *) myalloc (sizeof (OtherElem));
    other = &(other_elems->other_list[0]);
    other_elems->num_elems = 1;
  }
  else {
    other_elems = plyfile->other_elems;
    other_elems->other_list = (OtherElem *) realloc (other_elems->other_list,
                              sizeof (OtherElem) * (other_elems->num_elems+1));
    other = &(other_elems->other_list[other_elems->num_elems]);
    other_elems->num_elems++;
  }

  /* count of element instances in file */
  other->elem_count = elem_count;

  /* save name of element */
  other->elem_name = strdup (elem_name);

  /* create a list to hold all the current elements */
  other->other_data = (OtherData **)
                  malloc (sizeof (OtherData *) * other->elem_count);

  /* set up for getting elements */
  other->other_props = ply_get_other_properties (plyfile, elem_name,
                         offsetof(OtherData,other_props));

  /* grab all these elements */
  for (i = 0; i < other->elem_count; i++) {
    /* grab and element from the file */
    other->other_data[i] = (OtherData *) malloc (sizeof (OtherData));
    ply_get_element (plyfile, (void *) other->other_data[i]);
  }

  /* return pointer to the other elements data */
  return (other_elems);
}


void CPlyLoader::ply_get_element(PlyFile *plyfile, void *elem_ptr)
{
  if (plyfile->file_type == PLY_ASCII)
    ascii_get_element (plyfile, (char *) elem_ptr);
  else
    binary_get_element (plyfile, (char *) elem_ptr);
}


PlyOtherProp *CPlyLoader::ply_get_other_properties(
  PlyFile *plyfile,
  char *elem_name,
  int offset
)
{
  PlyElement *elem;
  PlyOtherProp *other;

  /* find information about the element */
  elem = find_element (plyfile, elem_name);
  if (elem == NULL) {
    fprintf (stderr, "ply_get_other_properties: Can't find element '%s'\n",
             elem_name);
    return (NULL);
  }

  other = get_other_properties (plyfile, elem, offset);
  return (other);
}


void CPlyLoader::close_ply(PlyFile *plyfile)
{
  fclose (plyfile->fp);
}

PlyFile *CPlyLoader::read_ply(FILE *fp)
{
  PlyFile *ply;
  int num_elems;
  char **elem_names;

  ply = ply_read (fp, &num_elems, &elem_names);

  return (ply);
}


PlyFile *CPlyLoader::ply_read(FILE *fp, int *nelems, char ***elem_names)
{
  int i,j;
  PlyFile *plyfile;
  int nwords;
  char **words;
  int found_format = 0;
  char **elist;
  PlyElement *elem;
  char *orig_line;

  /* check for NULL file pointer */
  if (fp == NULL)
    return (NULL);

  /* create record for this object */

  plyfile = (PlyFile *) myalloc (sizeof (PlyFile));
  plyfile->num_elem_types = 0;
  plyfile->comments = NULL;
  plyfile->num_comments = 0;
  plyfile->obj_info = NULL;
  plyfile->num_obj_info = 0;
  plyfile->fp = fp;
  plyfile->other_elems = NULL;
  plyfile->rule_list = NULL;

  /* read and parse the file's header */

  words = get_words (plyfile->fp, &nwords, &orig_line);
  if (!words || !equal_strings (words[0], "ply"))
    return (NULL);

  while (words) {

    /* parse words */

    if (equal_strings (words[0], "format")) {
      if (nwords != 3)
        return (NULL);
      if (equal_strings (words[1], "ascii"))
        plyfile->file_type = PLY_ASCII;
      else if (equal_strings (words[1], "binary_big_endian"))
        plyfile->file_type = PLY_BINARY_BE;
      else if (equal_strings (words[1], "binary_little_endian"))
        plyfile->file_type = PLY_BINARY_LE;
      else
        return (NULL);
      plyfile->version = atof (words[2]);
      found_format = 1;
    }
    else if (equal_strings (words[0], "element"))
      add_element (plyfile, words, nwords);
    else if (equal_strings (words[0], "property"))
      add_property (plyfile, words, nwords);
    else if (equal_strings (words[0], "comment"))
      add_comment (plyfile, orig_line);
    else if (equal_strings (words[0], "obj_info"))
      add_obj_info (plyfile, orig_line);
    else if (equal_strings (words[0], "end_header"))
      break;

    /* free up words space */
    free (words);

    words = get_words (plyfile->fp, &nwords, &orig_line);
  }

  /* create tags for each property of each element, to be used */
  /* later to say whether or not to store each property for the user */

  for (i = 0; i < plyfile->num_elem_types; i++) {
    elem = plyfile->elems[i];
    elem->store_prop = (char *) myalloc (sizeof (char) * elem->nprops);
    for (j = 0; j < elem->nprops; j++)
      elem->store_prop[j] = DONT_STORE_PROP;
    elem->other_offset = NO_OTHER_PROPS; /* no "other" props by default */
  }

  /* set return values about the elements */

  elist = (char **) myalloc (sizeof (char *) * plyfile->num_elem_types);
  for (i = 0; i < plyfile->num_elem_types; i++)
    elist[i] = strdup (plyfile->elems[i]->name);

  *elem_names = elist;
  *nelems = plyfile->num_elem_types;

  /* return a pointer to the file's information */

  return (plyfile);
}


void CPlyLoader::copy_property(PlyProperty *dest, PlyProperty *src)
{
  dest->name = strdup (src->name);
  dest->external_type = src->external_type;
  dest->internal_type = src->internal_type;
  dest->offset = src->offset;

  dest->is_list = src->is_list;
  dest->count_external = src->count_external;
  dest->count_internal = src->count_internal;
  dest->count_offset = src->count_offset;
}

PlyElement *CPlyLoader::find_element(PlyFile *plyfile, char *element)
{
  int i;

  for (i = 0; i < plyfile->num_elem_types; i++)
    if (equal_strings (element, plyfile->elems[i]->name))
      return (plyfile->elems[i]);

  return (NULL);
}


void CPlyLoader::add_element (PlyFile *plyfile, char **words, int nwords)
{
  PlyElement *elem;

  /* create the new element */
  elem = (PlyElement *) myalloc (sizeof (PlyElement));
  elem->name = strdup (words[1]);
  elem->num = atoi (words[2]);
  elem->nprops = 0;

  /* make room for new element in the object's list of elements */
  if (plyfile->num_elem_types == 0)
    plyfile->elems = (PlyElement **) myalloc (sizeof (PlyElement *));
  else
    plyfile->elems = (PlyElement **) realloc (plyfile->elems,
                     sizeof (PlyElement *) * (plyfile->num_elem_types + 1));

  /* add the new element to the object's list */
  plyfile->elems[plyfile->num_elem_types] = elem;
  plyfile->num_elem_types++;
}


void CPlyLoader::add_property (PlyFile *plyfile, char **words, int nwords)
{
  int prop_type;
//  int count_type;
  PlyProperty *prop;
  PlyElement *elem;

  /* create the new property */

  prop = (PlyProperty *) myalloc (sizeof (PlyProperty));

  if (equal_strings (words[1], "list")) {          /* list */
    prop->count_external = get_prop_type (words[2]);
    prop->external_type = get_prop_type (words[3]);
    prop->name = strdup (words[4]);
    prop->is_list = PLY_LIST;
  }
  else if (equal_strings (words[1], "string")) {   /* string */
    prop->count_external = Int8;
    prop->external_type = Int8;
    prop->name = strdup (words[2]);
    prop->is_list = PLY_STRING;
  }
  else {                                           /* scalar */
    prop->external_type = get_prop_type (words[1]);
    prop->name = strdup (words[2]);
    prop->is_list = PLY_SCALAR;
  }

  /* add this property to the list of properties of the current element */

  elem = plyfile->elems[plyfile->num_elem_types - 1];

  if (elem->nprops == 0)
    elem->props = (PlyProperty **) myalloc (sizeof (PlyProperty *));
  else
    elem->props = (PlyProperty **) realloc (elem->props,
                  sizeof (PlyProperty *) * (elem->nprops + 1));

  elem->props[elem->nprops] = prop;
  elem->nprops++;
}


void CPlyLoader::add_comment (PlyFile *plyfile, char *line)
{
  int i;

  /* skip over "comment" and leading spaces and tabs */
  i = 7;
  while (line[i] == ' ' || line[i] == '\t')
    i++;

  append_comment_ply (plyfile, &line[i]);
}


void CPlyLoader::append_comment_ply(PlyFile *ply, char *comment)
{
  /* (re)allocate space for new comment */
  if (ply->num_comments == 0)
    ply->comments = (char **) myalloc (sizeof (char *));
  else
    ply->comments = (char **) realloc (ply->comments,
		     sizeof (char *) * (ply->num_comments + 1));

  /* add comment to list */
  ply->comments[ply->num_comments] = strdup (comment);
  ply->num_comments++;
}


void CPlyLoader::add_obj_info (PlyFile *plyfile, char *line)
{
  int i;

  /* skip over "obj_info" and leading spaces and tabs */
  i = 8;
  while (line[i] == ' ' || line[i] == '\t')
    i++;

  append_obj_info_ply (plyfile, &line[i]);

}


void CPlyLoader::append_obj_info_ply(PlyFile *ply, char *obj_info)
{
  /* (re)allocate space for new info */
  if (ply->num_obj_info == 0)
    ply->obj_info = (char **) myalloc (sizeof (char *));
  else
    ply->obj_info = (char **) realloc (ply->obj_info,
		     sizeof (char *) * (ply->num_obj_info + 1));

  /* add info to list */
  ply->obj_info[ply->num_obj_info] = strdup (obj_info);
  ply->num_obj_info++;
}

int CPlyLoader::get_prop_type(char *type_name)
{
  int i;

  /* try to match the type name */
  for (i = StartType + 1; i < EndType; i++)
    if (equal_strings (type_name, type_names[i]))
      return (i);

  /* see if we can match an old type name */
  for (i = StartType + 1; i < EndType; i++)
    if (equal_strings (type_name, old_type_names[i]))
      return (i);

  /* if we get here, we didn't find the type */
  return (0);
}

