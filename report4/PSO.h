#pragma once

#include<vector>
#include <string>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

class PSO
{
private:
	int minNum = 0;
	double c1 = 0.5;
	double c2 = 0.5;
	double w = 0.5;
	std::vector<double> varMax, varMin;//変数の最小値・最大値
public:
	double resultSumValue;//評価関数の合計

	class Data//データ格納用クラス
	{
	public:
		std::vector<double> x;//座標
		std::vector<double> x_pbset;//座標
		std::vector<double> v;
		double functionValue;//与えられた関数の値
		double functionValuePbset;
		double result;

		Data(int _var_num) ://コンストラクタ
			x(std::vector<double>(_var_num)),
			v(std::vector<double>(_var_num)),
			x_pbset(std::vector<double>(_var_num))
		{
		}
	};

	std::vector<Data> data, prev_data;//操作前後で値を保持するために2個
	Data eliteData;
	PSO(int _max_genom_list, int _var_num, std::vector<double> _varMax, std::vector<double> _varMin);	//コンストラクタ
	bool selection();//選択
	void calc(bool enableDisplay);//評価関数の計算
	void setPosition(void);
private:
	void calcResult(bool enableSort);
	int random(int min, int max);
	double random(int min, double max);
	double random(double min, int max);
	double random(double min, double max);
	void displayValues();
public:
	~PSO();//デコンストラクタ
};