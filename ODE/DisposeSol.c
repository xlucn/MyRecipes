#include "NR.h"
/**
 * @brief free a struct type SODEsol
 * @param sol a struct type SODEsol
 */
void DisposeSODEsol(SODEsol sol)
{
	free(sol.t);
	for (int i = 0; i < sol.step; i++)
	{
		free(sol.y[i]);
	}
	free(sol.y);
}

/**
 * @brief free a struct type ODEsol
 * @param sol a struct type ODEsol
 */
void DisposeODEsol(ODEsol sol)
{
	free(sol.t);
	free(sol.y);
}
