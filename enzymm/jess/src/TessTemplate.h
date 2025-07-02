// ==================================================================
// TessTemplate.h
// Copyright (c) Jonathan Barker, 2002
// ==================================================================
// Declaration of TessTemplate creation.
// ==================================================================

#ifndef TESSTEMPLATE_H
#define TESSTEMPLATE_H

#include "Template.h"
#include <stdio.h>

// ==================================================================
// Creation of a TessTemplate object
// ==================================================================
// create(file,s)			Parse file & create template
// ==================================================================

extern Template *TessTemplate_create(FILE*,const char*);

// ==================================================================

#endif

