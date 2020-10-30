////////////////////////////////////////////////////////////////////////
// rtc_def.h
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
#ifndef RTC_RAYTRACE_DEF_H
#define RTC_RAYTRACE_DEF_H

namespace rtc { namespace cnst {// -------------------------------------
	// ���C�u�����������ŗ��p���Ă���萔�B/////////////////////
	enum plot_style {
		plot_xy,
		plot_xz,
		plot_yz,
	}; // �e���f���̃v���b�g���A�ǂ̎��ɉ����čs�����w�肷��B

	// �����萔 ////////////////////////////////////////////////
	const static double pi = 3.1415926535897932; // �~����
	const static double c  = 2.99792458e8;       // �����x[m/s]
	const static double e  = 1.60217759e-19;     // �d�C�f��
	const static double e0 = 8.85418782e-12;     // �U�d��
	const static double u0 = 4*pi*1e-7;          // ������
	const static double me = 9.1093922e-31;      // �d�q�̎���[kg]
	const static double mp = 1.6726217e-27;      // �z�q�̎���[kg]
	const static double k  = 1.3806505e-23;      // �{���c�}���萔

	const static double Rj = 71372.e3;           // �ؐ����a[m]

}}// namespace rtc::cnst; ----------------------------------------------

#endif//RTC_RAYTRACE_DEF_H
