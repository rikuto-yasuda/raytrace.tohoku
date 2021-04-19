// tracer.h: tracer �N���X�̃C���^�[�t�F�C�X
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRACER_H__05BCBAAA_04CB_4D2D_B090_5261F8726BEA__INCLUDED_)
#define AFX_TRACER_H__05BCBAAA_04CB_4D2D_B090_5261F8726BEA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "wave_param.h"

namespace rtc { // -----------------------------------------------------
	/* *********************************************************************
	+ ray
	+    �^����ꂽ���f���𗘗p���āA�g���[�X�v�Z�����s����B
	+   ���f���́Abasic_plasma_model��basic_magnetic_model��
	+   �h���N���X��n���Bray�́A�v�Z���ɂ��ꂼ��̃��f������
	+   �l���Q�Ƃ���B
	+    M �ɂ� basic_magnet_model����h���������ꃂ�f���N���X���A
	+   P�ɂ� basic_plasma_model����h�������v���Y�}���f���N���X��
	+   ���ꂼ��w�肵�Ȃ���΂Ȃ�Ȃ��B����ȊO�̃N���X���w�肷���
	+   �R���p�C�����ɃG���[�ƂȂ�B
	+
	+    �C���X�^���X�𐶐�������́Atake_a_step()���\�b�h���Ăяo�����Ƃ�
	+   �P�X�e�b�v���ʒu�A�g���x�N�g����i�߂邱�Ƃ��ł���B
	+   �P�X�e�b�v�ŕK���萔���Ԃ��i�ނƂ͌���Ȃ��̂Œ��ӂ��邱�ƁB
	+   �i�񂾎��Ԃ� take_a_step()�̖߂�l����m�邱�Ƃ��ł���B
	+
	***********************************************************************/
	class ray  
	{
	public:
		ray( const wave_parameter& wparam );
		virtual ~ray();

		// take_a_step()�͈ʒu�Ɣg���x�N�g����1�X�e�b�v�i�߂郁�\�b�h�ł���B
		// r �ɂ͌��݂̔g���ʒu���Ak �ɂ͌��݂̔g���x�N�g�����w�肷��B
		// �Ăяo��������ɁAr��k���V�����l�ɍX�V�����B
		// ��{�I�ɃR���X�g���N�^�ɓn���� timeStep���Ԍo�ߌ��r,k��
		// ���ʂƂ��ĕԂ����B�������A�������ɂ߂ĕ��G�ȓ���������\��������ꏊ�ł�
		// ����ɍׂ������Ԍo�ߌ�̌��ʂ��Ԃ���邱�Ƃ�����B
		// ���ۂɌv�Z���ꂽ�P�X�e�b�v�Ԃ̎��ԊԊu�́A�߂�l�ŕԂ����B
		//
		// ��ڂ̈����́A�w�肵��vector_pair��r��k�x�N�g���̍���(dr, dk)��
		// �i�[����Đ��䂪������B�ȗ��B
		inline double take_a_step();

		// �ʒu���̎擾 ////////////////////////
		// ���݂̔g���ʒu���ASM���W�n�ŕԂ��B
		const vector& getR() const
		{ return m_rk.first; }

		// ���݂̔g���x�N�g�����ASM���W�n�ŕԂ��B
		const vector& getK() const
		{ return m_rk.second; }

		// ���O�̈ʒu�Ƃ̍����ASM���W�n�ŕԂ��B
		const vector& getDeltaR() const
		{ return m_drk.first; }

		// ���O�̔g���x�N�g���Ƃ̍����ASM���W�n�ŕԂ��B
		const vector& getDeltaK() const
		{ return m_drk.second; }


		// �����̔g���x�N�g���𐶐����邽�߂̕⏕�֐��B
		ray* initialize(
			const vector& r,
			const vector& k
		);
		ray* initialize(
			double rx, double ry, double rz,
			double kx, double ky, double kz
		);
		ray* initialize(
			const vector& r,
			double pitch, double round
		);
		ray* initialize(
			double rx, double ry, double rz,
			double pitch, double round
		);

		// ���N���X���擾����B
		wave_parameter&       getWaveParam()       { return m_wave; }
		const wave_parameter& getWaveParam() const { return m_wave; }

		// ���O���b�Z�[�W�̎擾
		const std::ostringstream& getLog() const
		{ return m_log; }
		
	private:
		// �v�Z�⏕�֐� //
		// �x�[�X�ƂȂ鎮�́Areference�f�B���N�g���ȉ��̉摜�t�@�C�������邱�ƁB

		// intermediate�\���̂ɂ́A�v�Z�r���ł悭���p����l��
		// ���炩���ߌv�Z���Ċi�[���Ă����B
		// 
		// m_im �� makeInitialVector()�����ŏ���������A
		// take_a_step()�ŗ��p���A�����Ԃ��O��
		// update_intermediate()���Ăяo���A���̃X�e�b�v�p�ɍX�V����B
		//
		// �v�Z����update_intermediate()�֐��������݂邱�ƁB
		//
		struct intermediate
		{
			vector B;
			double
				Bk,
				B2,
				k2,
				w,
				w2,
				wp2,
				wc2,
				s,c,
				s2,c2,
				X2,Y2, X4, Y4,
				iX2, iX22,
				numerator,
				denominator,
				root
			;

			intermediate& operator = ( const intermediate& r );
			bool operator ==( const intermediate& r ) const;

		} m_im;

		void update_intermediate(
			intermediate& i,
			const vector& r,
			const vector& k
		) const ;

		// 1step�̎��ԗʂ��v�Z���A�ړ��ɂ����鎞�Ԃ�
		// �ړ���̓_�ɂ����� intermediate ��Ԃ��B
		double calc_dt(
			const vector_pair& rk,
			const vector_pair& drk,
			intermediate&      out_i
		) const;

		// �ꎟ������ ///////////////////////////
		double calc_dGdw(
			const intermediate& i,
			const vector&       r,
			const vector&       k
		) const; // dG/dw���v�Z���A���ʂ�Ԃ��B

		vector calc_dGdr(
			const intermediate& i,
			const vector&       r,
			const vector&       k
		) const; // dG/dr���v�Z���A���ʂ�Ԃ��B

		vector calc_dGdk(
			const intermediate& i,
			const vector&       r,
			const vector&       k
		) const; // dG/dk���v�Z���A���ʂ�Ԃ��B

		// �񎟔����� ///////////////////////////
		double numerator_G(
			const intermediate& i,
			const vector&       r
		)   const; // �֐�G�̑�Q���̕��q�̒l��Ԃ��B
		double denominator_G(
			const intermediate& i,
			const vector&       r
		) const; // �֐�G�̑�Q���̕���̒l��Ԃ��B

		int reflectioncheck(
			const vector&       r,
			const vector_pair& drk
		) const;

		double reflect_dt(
			const vector_pair rk,
			const vector_pair drk,
			intermediate&      out_i
		) const;

		vector reflect_n(
			const vector_pair rk,
			intermediate&      out_i
		) const;


		// �G���[�`�F�b�N ///////////////////////
		void checkState(
			const intermediate& i,
			const vector_pair& drk,
			const vector& r,
			const vector& k
		) const;

		// ���O���b�Z�[�W�̑}�� /////////////////
		const char* log( const char* msg ) const
		{ m_log << msg << "\n"; return msg; }

	private:
		// �ʒu /////////////////////////////////
		vector_pair m_rk, m_drk;

		// �p�����[�^ ///////////////////////////
		wave_parameter  m_wave;
		
		// dt�v�Z�ɗ��p����p�����[�^
		// mkInitialVector()�ōŏ��l�ɏ���������A
		// calc_dt()���Ăяo����邽�тɍX�V����B
		mutable double m_dt_before; 
		
		// �g���[�X�E���O ///////////////////////
		mutable std::ostringstream m_log;
	};
}// namespace rtc; -----------------------------------------------------

#include "tracer.inl"

#endif // !defined(AFX_TRACER_H__05BCBAAA_04CB_4D2D_B090_5261F8726BEA__INCLUDED_)
