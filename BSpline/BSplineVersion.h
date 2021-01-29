#ifndef _BSPLINE_VERSION_H_
#define _BSPLINE_VERSION_H_

/*
 * Provide a default revision which will be exact for tagged release builds,
 * and otherwise is generic unless the automatic revision info is provided
 * by a build.  The generic form follows the automatic form of v1.6-count,
 * except the count is replaced with x.  The other common form would be the
 * semantic version prerelease, like v1.7-alpha, except that seems confusing
 * when compared with the git revisions.  So prerelease tags can still be
 * used, but they will have to be explicitly defined in this file,
 * committed, and tagged that way.  The following commits would use the form
 * v1.7-alpha-x.  Ideally, the auto revisions are only seen by developers
 * and never appear in releases.
 */
#ifdef BSPLINE_AUTO_REVISION
// We have to use two levels of macros to first expand BSPLINE_AUTO_REVISION
// before stringifying it.
#define bspline_xstr(s) bspline_str(s)
#define bspline_str(s) #s
#define BSPLINE_VERSION  bspline_xstr(BSPLINE_AUTO_REVISION)
#else
#define BSPLINE_VERSION  "v1.6-x"
#endif

/*
 * The repo URL has been hardcoded to the NCAR organization, since it moved
 * from the ncareol organization.
 */
#define BSPLINE_URL      "https://github.com/NCAR/bspline"

#endif // _BSPLINE_VERSION_H_
