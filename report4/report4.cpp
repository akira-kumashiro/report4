#include "stdafx.h"
#include "PSO.h"

#define MAX_GENERATION 30000
#define MAX_GENOM_LIST 50
#define MAX_WEIGHT 60
#define VAR_NUM 2
//min:(-0.54719,-1.54719)=-1.9133

int main()
{
	double varMax[] = { 4,4 };
	double varMin[] = { -1.5,-3 };

	//配列をstd::vectorへ変換
	std::vector<double> vMax(varMax, std::end(varMax));
	std::vector<double> vMin(varMin, std::end(varMin));
	PSO pso(MAX_GENOM_LIST, VAR_NUM, vMax, vMin);//遺伝的アルゴリズム諸関数をまとめたクラスの宣言

	for (int i = 0; i <= MAX_GENERATION; i++)//メインのループ
	{
		bool change = pso.selection();//選択

		if (i % (MAX_GENERATION / 10) == 0 || change)
		{
			std::cout << "i=" << std::to_string(i) << std::endl;
			pso.calc(true);//評価関数の計算
		}
		else
		{
			pso.calc(false);//評価関数の計算
		}
	}

	while (1)
	{
		if (_kbhit() && _getch() == 27)
			break;
	}
	return 0;
}
