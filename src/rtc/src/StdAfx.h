// stdafx.h : �W���̃V�X�e�� �C���N���[�h �t�@�C���A
//            �܂��͎Q�Ɖ񐔂������A�����܂�ύX����Ȃ�
//            �v���W�F�N�g��p�̃C���N���[�h �t�@�C�����L�q���܂��B
//

#if !defined(AFX_STDAFX_H__89C5ABD0_D5FE_4AB7_8153_71EF8F228A65__INCLUDED_)
#define AFX_STDAFX_H__89C5ABD0_D5FE_4AB7_8153_71EF8F228A65__INCLUDED_

#if defined _WIN32 && _MSC_VER > 1000
#pragma once

#define WIN32_LEAN_AND_MEAN		// Windows �w�b�_�[����w�ǎg�p����Ȃ��X�^�b�t�����O���܂�


// TODO: �v���O�����ŕK�v�ȃw�b�_�[�Q�Ƃ�ǉ����Ă��������B

//��MFC�Ń��������[�N�񍐂���ɂ́A�ȉ��̃R�[�h��ǉ�����
//main()�̐擪�Ɂu_CrtSetDbgFlag( _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG) | _CRTDBG_LEAK_CHECK_DF);�v�������B
#ifdef _DEBUG
#include <cstdlib>
#include <new>
#include <memory>

#include <crtdbg.h>
#define _CRTDBG_MAP_ALLOC

#define new  ::new(_NORMAL_BLOCK, __FILE__, __LINE__)
#endif
#endif // _MSC_VER > 1000

// �W�����C�u���� ///////////////////////////////
#include <cstdlib>
#include <cfloat>

#include <vector>
#include <list>
#include <set>

#if defined __INTEL_COMPILER
#include "../iccmath.hpp"
#else
#include <cmath>
#include <complex>
#endif

// boost library ////////////////////////////////
#include <boost/format.hpp>
#include <boost/math/quaternion.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/thread.hpp>

// MSVC6�ł͎g���Ȃ�
#if (defined (_MSC_VER) && _MSC_VER > 1200) || !defined _MSC_VER
#	include <boost/numeric/ublas/lu.hpp>
#endif

// namespace ////////////////////////////////////
#if defined (_MSC_VER) && _MSC_VER < 1300
namespace std {
	using ::memset;
	using ::memcpy;
}
#endif

// libraytrace //////////////////////////////////
#include "../raytrace.h"

// STLport�̃`�F�b�N ///////////////////////////////////////////////////
#ifndef NDEBUG
#	if !defined _STLP_USE_DYNMIC_LIB || !defined _STLP_DEBUG
#	pragma message("DEBUG�ł� _STLP_USE_DYNMIC_LIB �� _STLP_DEBUG ����`����ĂȂ���H")
#	endif
#endif

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ �͑O�s�̒��O�ɒǉ��̐錾��}�����܂��B

#endif // !defined(AFX_STDAFX_H__89C5ABD0_D5FE_4AB7_8153_71EF8F228A65__INCLUDED_)
