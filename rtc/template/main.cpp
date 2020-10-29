////////////////////////////////////////////////////////////////////////
// RTC template main.cpp
//  This program has been written in C++.
//  Copyright (C) 2006 Miyamoto Luisch.

// include ���� ����ȊO�������Ă͂����܂���B
// �ǉ��� include �́A���ׂ�stdlib.h�ɏ����Ă��������B
#include "stdlib.h"

// NOTE:
// ����́ARTC��p�����R�[�h�������Ƃ��́u�ЂȌ`�v�ƂȂ�悤��
// �K�v�Œ���̃R�[�h���L�q���ꂽ�e���v���[�g�ł��B
// RTC���g�p����ꍇ�Atemplate �ɂ���t�@�C�������̂܂܃R�s�[����
// �������ɕύX�������Ă����ƕ֗��ł��B
//
// �R���p�C���Ő����������s�\�t�@�C���̖��O��
//  common.mak ���� LIBNAME = ���Œ�`����Ă��܂��B
//
// makefile�́A�R���p�C���ɉ�����
//  GNU gcc              : makefile.gcc
//  Microsoft Visual C++ : makefile.vc
//  Intel Compiler       : makefile.icc
// ���g�p���Ă��������B�R�}���h���C������R���p�C������Ƃ��́A
//
//  $ make -f makefile.gcc
//
// �̗l�Ɏg���܂��BUNIX��Ŏg�p����ꍇ�́A�Ⴆ��
//
//  $ ln -s makefile.gcc makefile
//
// �̗l�� makefile����V���{���b�N�����N�𒣂�ƕ֗��ł��B
//
// �����炵���\�[�X�t�@�C�����������ꍇ�́Acommon.mak����
// SRC = ���ɒǉ����Ă��������B
// �܂��A��Ƀ\�[�X�R�[�h�̐擪�� stdlib.h ���C���N���[�h����悤
// �C�����Ă��������B��������Ȃ��ꍇ�A�ꕔ�̃R���p�C����
// �R���p�C�������s���܂��B


// main �֐� -------------------------------------
int main()
{
	// RTC�ł́A�g���[�X���̃G���[�����ɂ͗�O������p���܂��B
	// ���ׂẴG���[���m���ɕ⑫���邽�߁Amain�̐擪�� try�߂������܂��B
	try {
		// rtc::cosmos �̂݁A�O���[�o���ϐ��ɒu���Ă����܂��܂���B
		// ���A�\�z���ɃG���[�����������ꍇ�uSegmentation fault�v���������܂��B
		// try ���̒��ɋL�q���鎖�ŁA�����������x�h�������ł��܂��B
		rtc::cosmos c(1980,5,6, 0,0,0);
		
		// ����ȍ~�Artc::getCosmos() �� ���O�ɍ쐬����
		// cosomos �I�u�W�F�N�g��������悤�ɂȂ�܂��B
		assert( &c == &rtc::getCosmos() );
		
		// ���f���́A����̕����܂߁A�O���[�o���ϐ��ɒu���Ȃ��ł��������B
		// �O���[�o���ō\�z�����ꍇ�A�v���I�G���[���������邱�Ƃ�����܂��B
		// �Ȃ��A���f���̓��C�g���[�X�����s����ԁA�j�󂵂Ȃ��ł��������B
		// planet�̓��f�����R�s�[�����A�^����ꂽ�C���X�^���X�𒼐ڎg�p���܂��B
		rtc::model::magnet::simple m;
		rtc::model::plasma::simple p;
		
		// �f���́A����̕����܂߁A�O���[�o���ϐ��ɒu���Ȃ��ł��������B
		// �O���[�o���ō\�z�����ꍇ�A�v���I�G���[���������邱�Ƃ�����܂��B
		rtc::planet::earth e(m,p);
		
		// ���������f�����F���ɓo�^���܂��B
		c.registerPlanet(e);
		
		// ���C�g���[�X�̏����ʒu�͎��̃R�[�h�Ŏ擾�o���܂��B
		rtc::vector start_ptr = e.getLocation( 70.0/*ILAT*/, 0.0/*MLT*/, 3e7/*Altitude[m]*/ );
		
		// ���C�g���[�X���邽�߂ɂ́u���v�����K�v������܂��B
		// RTC�ł�cosmos�� ray �I�u�W�F�N�g�𐶐����A���̃|�C���^��Ԃ������ł��܂��B
		// ���̃R�[�h�ŁA400[kHz]��R-X���[�h�g�����쐬���邱�Ƃ��ł��܂��B
		rtc::ray* r = c.createRay( rtc::wave_parameter::RX_MODE, 400e3 );
		
		// ���̃R�[�h�Ő������ꂽ���̌����������x�N�g���ik�x�N�g���j��ݒ肷�邱�Ƃ��ł��܂��B
		r->initialize( start_ptr, 90.0/*pitch*/, 180.0/*round*/ );
		
		// 5000�X�e�b�v�����g���[�X���A���̌��ʂ�W���o�͂ɏo���܂��B
		for( int n = 0; n < 5000; ++n )
		{
			// ray::take_a_step()�ŁA�K�v�Ȉ�X�e�b�v�i�߂鎖���ł��܂��B
			// ray::take_a_step()�́A��X�e�b�v�ԂɌo�߂���������Ԃ��܂��B
			r->take_a_step();
			
			// ray �̌��݈ʒu�� ray::getR()�Ŏ擾���鎖���ł��܂��B
			const rtc::vector& pos = r->getR();
			std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
		}
		
		// �����������́A���̃R�[�h�Ō�n�����邱�Ƃ𐄏����܂��B
		// �K�{�ł͂���܂���B
		c.eraseRay(r);
	}
	catch( std::exception& e )
	{
		// RTC�������瓊����ꂽ��O�G���[�͂����ŃL���b�`����
		// �G���[�̓��e��\�����ďI�����܂��B
		std::cerr << e.what() << std::endl;
	}

	return 0;
}
