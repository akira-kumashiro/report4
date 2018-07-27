#include "stdafx.h"
#include "PSO.h"

PSO::PSO(int _max_genom_list, int _var_num, std::vector<double> _varMax, std::vector<double> _varMin) :
	data(std::vector<Data>(_max_genom_list, _var_num)),//dataの初期化
	eliteData(_var_num)
{
	//もらった変数をクラス内変数に格納
	varMax = _varMax;
	varMin = _varMin;

	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			data[i].x[j] = random(varMin[j], varMax[j]);//遺伝子の初期設定
		}
	}
	prev_data = data;
	calcResult(true);

	displayValues();
}

bool PSO::selection()
{
	int max_num = 0;//最も評価の良い個体の番号
	bool ret = false;

	calc(false);

	for (int i = 0; i < data.size(); i++)//ルーレット選択用に評価関数の合計と一番評価の良い番号を取得
	{
		if (data[i].result > data[max_num].result)
			max_num = i;
	}

	eliteData = data[max_num];//最も評価の良い個体を保持
	if (prev_data[max_num].functionValue - eliteData.functionValue != 0)//最も評価の良い個体の変化の監視(デバッグ用)
		ret = true;

	setPosition();
	return ret;
}

void PSO::calc(bool enableDisplay)
{
	calcResult(true);
	for (int i = 0; i < data.size(); i++)//評価関数が最小の奴と最大のやつを検索
	{
		if (data[i].result < data[minNum].result)
			minNum = i;
	}
	//評価関数が最もいいやつを保存
	data[minNum] = eliteData;

	calcResult(true);

	if (enableDisplay)
		displayValues();
}

void PSO::setPosition(void)
{
	prev_data = data;
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].v.size(); j++)//位置・速度設定
		{
			data[i].v[j] = w * data[i].v[j] + c1 * random(0.0, 1.0)*(data[i].x_pbset[j] - data[i].x[j]) + c2 * random(0.0, 1.0)*(eliteData.x[j] - data[i].x[j]);
			data[i].x[j] += data[i].v[j];
		}
	}
}

void PSO::calcResult(bool enableSort)
{
	int maxNum = 0;
	double seg;

	for (int i = 0; i < data.size(); i++)
	{
		data[i].functionValue = std::sin(data[i].x[0] + data[i].x[1]) + std::pow((data[i].x[0] - data[i].x[1]), 2.0) - 1.5*data[i].x[0] + 2.5*data[i].x[1] + 1;//与えられた関数
		if (data[i].functionValue < data[i].functionValuePbset)
		{
			data[i].functionValuePbset = data[i].functionValue;
			data[i].x_pbset = data[i].x;

		}
	}
	seg = data[maxNum].functionValue;//評価関数の切片を与えられた関数が最も大きいやつにセット
	resultSumValue = 0;

	for (int i = 0; i < data.size(); i++)
	{
		bool flag = true;
		double coefficient = 0.01;//評価関数用の定数
		for (int j = 0; j < data[i].x.size(); j++)
		{
			if (data[i].x[j] > varMax[j] || data[i].x[j] < varMin[j])//座標が場外にいるやつの処理
				flag = false;
		}
		data[i].result = std::pow((data[i].functionValue - seg), 2.0);//与えられた関数の値から切片で設定した値を引いて2乗する→与えられた関数の値が小さいやつが強くなる

		if (!flag)//場外に出たやつの処理
			data[i].result *= coefficient;
		resultSumValue += data[i].result;
	}
	if (enableSort)
		std::sort(data.begin(), data.end(), [](const Data& x, const Data& y) { return x.functionValue > y.functionValue; });
}

int PSO::random(int min, int max)
{
	//乱数の設定
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_int_distribution<int> distribution(min, max);
	return distribution(engine);
}

double PSO::random(int min, double max)
{
	return random((double)min, max);
}

double PSO::random(double min, int max)
{
	return random(min, (double)max);
}

double PSO::random(double min, double max)
{
	std::random_device rnd;
	std::mt19937 engine(rnd());
	std::uniform_real_distribution<double> distribution(min, max);
	return distribution(engine);
}

void PSO::displayValues()
{
	for (int i = 0; i < data.size(); i++)
	{
		for (int j = 0; j < data[i].x.size(); j++)
		{
			printf_s("%10.7lf,", data[i].x[j]);//デバッグ用
		}
		printf_s(" \t f(x,y)=%10.7lf\t Result=%10.7lf\n", data[i].functionValue, data[i].result);
	}
}

PSO::~PSO()
{

}