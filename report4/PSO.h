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
	int max_genom_list;//�̐�
	int var_num;//�i���̌�
	int minNum = 0, maxNum = 0;
	double c1 = 0.5;
	double c2 = 0.5;
	double w = 0.5;
	std::vector<double> varMax, varMin;//�ϐ��̍ŏ��l�E�ő�l
public:
	double resultSumValue;//�]���֐��̍��v

	class Data//�f�[�^�i�[�p�N���X
	{
	public:
		std::vector<double> x;//���W
		std::vector<double> x_pbset;//���W
		std::vector<double> v;
		double functionValue;//�^����ꂽ�֐��̒l
		double functionValuePbset;
		double result;

		Data(int _var_num)://�R���X�g���N�^
			x(std::vector<double>(_var_num)),
			v(std::vector<double>(_var_num)),
			x_pbset(std::vector<double>(_var_num))
		{
		}
	};

	std::vector<Data> data, prev_data;//����O��Œl��ێ����邽�߂�2��
	Data eliteData;
	PSO(int _max_genom_list, int _var_num, std::vector<double> _varMax, std::vector<double> _varMin);	//�R���X�g���N�^
	bool selection();//�I��

	//bool blxAlphaCrossover();
	//bool mutation();//�ˑR�ψ�
	bool calc(bool enableDisplay);//�]���֐��̌v�Z
	void setPosition(void);
private:
	bool calcResult(bool enableSort);
	int random(int min, int max);
	double random(int min, double max);
	double random(double min, int max);
	double random(double min, double max);
	bool displayValues();
	; public:
		~PSO();//�f�R���X�g���N�^
};