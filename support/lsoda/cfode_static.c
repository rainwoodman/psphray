#include "common.h"
#include "cfode_static.inc"

void cfode_static (int meth)
{
#ifdef CFODE_STATIC
	if (meth == 1) {
		_C(elco) = elco1;
		_C(tesco) = tesco1;

	} else {
		_C(elco) = elco2;
		_C(tesco) = tesco2;
	}
#endif
}
