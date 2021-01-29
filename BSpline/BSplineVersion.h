#ifndef _BSPLINE_VERSION_H_
#define _BSPLINE_VERSION_H_

/*
 * Provide a default revision which will be exact for tagged release builds,
 * and otherwise is generic unless the automatic revision info is provided
 * by a build.
 */
#define BSPLINE_VERSION  "v1.7-devel"

/*
 * The repo URL has been hardcoded to the NCAR organization, since it moved
 * from the ncareol organization.
 */
#define BSPLINE_URL      "https://github.com/NCAR/bspline"

#ifdef BSPLINE_AUTO_REVISION
#include "bspline-auto-revision.h"
#endif

#endif // _BSPLINE_VERSION_H_
