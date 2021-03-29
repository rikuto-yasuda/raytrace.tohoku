////////////////////////////////////////////////////////////////////////
// libraytrace raytrace.cpp
//  This program has been written in C++.
//  Copyright (C) 2006 Miyamoto Luisch.
#include "testing.h"
#include "raytrace.h"

// -------------------------------------------------------------------
// raytrace::raytrace()
// ���l�ɏ]���ď����l���쐬���Aoperator()�ɐ����n����
// �g���[�X�����s����Brtc::ray�̓R���X�g���N�^�Ő�������B
//
raytrace::raytrace(
	const testing_env*   env,
	const double       round
) : m_env( env  ),
    m_ray( NULL ),
	m_round(round),
	m_progress( 0.0 ),
    m_state("init")
{
	// �����\�z���邽�߂ɁArtc::cosmos::createRay()���Ăяo��
	// �F���Ɍ��̃C���X�^���X�𐶐�������B
	// ����������̃p�����[�^�́A���ׂĈ����Ŏw�肷��B
	rtc::ray* ray = rtc::getCosmos().createRay(
		env->mode,                   // LO_MODE��RX_MODE���w�肷��B
		(2*rtc::cnst::pi)*env->freq, // �g���̊p���g����[Hz]�Ŏw�肷��B
		env->precision,              // �g��step�O��ł̋��e�����w�肷��B�ȗ��B
		env->step_length,            // 1step���ɐi�߂�ő咷��[m]�Ŏw�肷��B�ȗ��B
		env->time_range.max,         // 1step���ɐi�߂鎞�Ԃ̍ő�l���w�肷��B�ȗ��B
		env->time_range.min          // 1step���ɐi�߂鎞�Ԃ̍ŏ��l���w�肷��B�ȗ��B
	);
	m_ray = ray;
	
	// ���C�p�X�̌����m��
	m_raypath.reserve( env->ray_segment );
///////////////////////////////////////////////////////step���m�F�p
	m_rayvariation.reserve( env->ray_segment );
	
	if( env->is_verbose )
	{
		// ���E���h�p�̕������o��
		m_output << "## round : " << rtc::rad2deg(round) << " ##\n";
	}
}


// -------------------------------------------------------------------
// raytrace::~raytrace 
// raytrace�N���X�̃f�X�g���N�^�B
// �������ꂽ ray �ւ̃|�C���^����������B
//
raytrace::~raytrace()
{
	// �Ō�ɁA�쐬�������͏�������B
	rtc::getCosmos().eraseRay( m_ray );
}


// -------------------------------------------------------------------
// raytrace::operator ()()
// ������������Ƀ��C�����[�v���Ăяo���B
// �I����Ƀo�b�N�g���[�X�̃`�F�b�N�����A
// �K�v�Ȃ�Ώ��������s������Ƀo�b�N�g���[�X����B
//
void raytrace::operator ()()
{
	// initialize()�͗�O�G���[�𓊂���\��������B
	// mainloop()��������O�𓊂���\�������邩��A
	// �������� catch �����`����B
	try
	{
		// �����������̈ʒu�x�N�g���A�g���x�N�g�������������Ȃ���΂Ȃ�Ȃ��B
		// �����̃x�N�g���y�A�́A�ʒu�ƌ������w�肵��
		// rtc::ray::initialize()���Ăяo�����Ƃ�
		// �ȒP�ɏ��������邱�Ƃ��ł���B
		//
		// initialize()�̂U�����łł́A�������W�n��
		// �����ʒu�Ɣg���x�N�g���̌������w�肷��Binitialize()��
		// �����f������K�؂Ȓl�ɕϊ����A�x�N�g���y�A���������B
		//
		// �R�����łł́A�s�b�`�p�ƃ��E���h�p���w�肷�邱�Ƃ�
		// ���͐��ɑ΂���C�ӂ̊p�x�����邱�Ƃ��ł���B
		m_ray->initialize(
			m_env->source_x,
			m_env->source_y,
			m_env->source_z,
			m_env->pitch_angle,
			m_round
		);
		m_state = "run";
		mainloop();
		
		// �o�b�N�g���[�X ----------------------------
		if( m_env->is_back_trace )
		{
			m_state = "back";

			m_output << "\n";
			if( m_env->is_verbose ) {
				m_output << "# back trace\n";
			}

			//
			// raytrace_proc()�̏I�_�ʒu����Ark�x�N�g�����t�����Ɍ�����
			// �������ŋt�����Ƀg���[�X����B
			//
			// �K�� makeInitialVector() ���Ăяo���ď��������邱�ƁB
			// �������Ȃ������ꍇ�̃g���[�X���ʂ͖���`�ƂȂ�B
			//
			const rtc::vector_pair rk = m_raypath.back().second;
			m_ray->initialize(
				 rk.first, // r
				-rk.second // k
			);

			m_raypath.clear();
			mainloop();
		}
		m_state = "comp.";
	}
	
	// �G���[�͎��s���̏����G���[�̑��ɁA
	// �������Ȃǂ̗��R�ł���ȏ�g���[�X�ł��Ȃ�����
	// �����Ă��邱�Ƃ�����B
	// ��{�I�ɃG���[���b�Z�[�W��\�����ďI�����邩
	// ���̌����ɏ������ڂ��΂悢�B

	catch ( std::range_error& e)
	{
		if( m_env->is_verbose ) {
			m_output << "## error ## " << e.what() << std::endl;
		}
		m_state = e.what();
	}

	catch ( std::exception& e )
	{
		if( m_env->is_verbose ) {
			m_output << "## error2 ## " << e.what() << std::endl;
		}
		m_state = e.what();	
	}

	catch( ... )
	{
		m_state = "error";
		throw;
	}
}
// -------------------------------------------------------------------
// raytrace::mainloop()
// �w�肳�ꂽ���A�����A�����ʒu����g���[�X�����s����B
// ���l�Ɏw�肳�ꂽ�I�_�ɂ��ǂ蒅�����Ƃ��A�����Ԃ��B
// ��{�I�Ƀ��[�L���O�X���b�h�œ��삳����̂ŁA���ӂ��邱�ƁB
//
void raytrace::mainloop()
{
	// ���̈ʒu���v���b�g���鋗���Ԋu�B
	const double ray_segment
		= m_env->ray_segment != 0.0 ?
		m_env->ray_length / m_env->ray_segment : -1.0;

	// ���s�����_�܂ł̃��C�p�X�̒������i�[����B
	double ray_length      = 0.0;

	// �g����������̌o�ߎ���[s]���i�[����B
	double t = 0.0;

	// rtc::ray::take_a_step()���Ăяo�����Ƃ�
	// 1step���i�܂��邱�Ƃ��ł���B
	// take_a_step()�́A�X�e�b�v�Ԃ̌o�ߎ��Ԃ�Ԃ��B
	// �g���̌��݈ʒu�́Aray::getR()�̖߂�l���瓾����B
	// ray::getDeltaR()����O��̍����x�N�g����������̂ŁA
	// ����𑫂����킹�A���H���� m_env->ray_length �𒴂����Ƃ��ɏI������B
	// �܂��A���H���� m_env->ray_length/m_env->ray_segment (== ray_segment)
	// �ɒB���閈�Ƀ��|�[�g�\�����s���B
	
	// rtc::ray::take_a_step()�͗�O�𓊂��Ă��邱�Ƃ�����B
	// �����߂܂��Ȃ��Ɓu�s���ȏ����������v����̂ŁA
	// try{} catch(...)�Ŋm���ɕ߂܂��A�G���[�������s���B
	// ���̏����� operator () �ɋL�q���Ă���̂ŎQ�Ƃ̂��ƁB
	
	// ���߂ɊJ�n�n�_��\������B
	if( m_env->is_plot_startptr )
	{
		report_progress( .0 );
		print_location(raypath_element(
			0.0,
			rtc::vector_pair( m_ray->getR(), m_ray->getK() )
		));
////////////////////////////////////////////////////////////////step�Ԋu�m�F�p
		print_variation(raypath_element(
			0.0,
			rtc::vector_pair( m_ray->getR(), m_ray->getK() )
		));
///////////////////////////////////////////////////////////////
	}
	
	for( unsigned loop = 0;
		 loop < m_env->step_count;
		 ++loop
	){
		const double dt = m_ray->take_a_step(); // �g����1step�i�߂�b
		t += dt;
		
		const double dr = rtc::norm_2( m_ray->getDeltaR() );
		ray_length += dr;

		// �ȉ��A�����i�񂾋�������I���_�𓱏o�B
		if( ray_length < m_env->ray_length )
		{
			
			// �v���O���X�E�o�[��\��
			report_progress( std::max(
				ray_length / m_env->ray_length,
				static_cast<double>(loop) / m_env->step_count
			));
			
			// ���ʂ����߂�
			m_raypath.push_back(raypath_element(
				t,
				rtc::vector_pair( m_ray->getR(), m_ray->getK() )
			));
//////////////////////////////////////////////////////////////////////step�Ԋu�m�F�p
			m_rayvariation.push_back(raypath_element(
				dt,
				rtc::vector_pair( m_ray->getDeltaR(), m_ray->getDeltaK() )
			));
///////////////////////////////////////////////////////////////////////
		}
		else /*( ray_length >= m_env->ray_length )*/
		{
			//�ŏI�_��ray_length���傤�ǂɎ��߂�B
			rtc::vector r = m_ray->getR()-m_ray->getDeltaR();
			rtc::vector k = m_ray->getK()-m_ray->getDeltaK();
			
			const double lest = m_env->ray_length - (ray_length-dr);
			const double factor = lest/dr;
			r += factor * m_ray->getDeltaR();
			k += factor * m_ray->getDeltaK();
			
			m_raypath.push_back(raypath_element(
				t - dt*(1.0-factor),
				rtc::vector_pair( r, k )
			));
///////////////////////////////////////////////////////////////////////step�Ԋu�m�F�p
			m_rayvariation.push_back(raypath_element(
				dt * factor,
				rtc::vector_pair( factor * m_ray->getDeltaR(), factor * m_ray->getDeltaK() )
			));
///////////////////////////////////////////////////////////////////////
			break;
		}
	}

	// m_output�ɐ������ďo��
	const double n = m_raypath.size();
	if( n < m_env->ray_segment )
	{
		int i = 0;
		for( i = 0; i < n; ++i )
		{
			print_location( m_raypath[i] );
/////////////////////////////////////////////////////////step�Ԋu�m�F�p
			print_variation( m_rayvariation[i] );
		}
		for( ; i < m_env->ray_segment; ++i )
		{
			print_location( m_raypath[ static_cast<int>(n-1) ] );
///////////////////////////////////////////////////////////////////////
			print_variation( m_rayvariation[ static_cast<int>(n-1) ] );
		}
	}
	
	else for( int i = 0; i < m_env->ray_segment; ++i )
	{
		print_location(
			m_raypath[ static_cast<int>(i * n/m_env->ray_segment) ]
		);
		print_variation(
			m_rayvariation[ static_cast<int>(i * n/m_env->ray_segment) ]
		);
	}
	
	// �I���n�_���m���ɕ\��������B
	report_progress( 1.0 );
	print_location( m_raypath.back() );
	print_variation( m_rayvariation.back() );
}


// -------------------------------------------------------------------
// raytrace::report_progress()
// �g���[�X���̃v���O���X�E�o�[��\������
//
void raytrace::report_progress( const double percent )
{
	//�i���󋵂��
	m_progress = percent;

	if( !m_env->is_parallel && 0 == (static_cast<int>(100*percent) % 3) )
	{
		std::clog << "[";

		int n = 0;
		for( ; n < 24*percent; ++n ) std::clog << "=";
		for( ; n < 24;         ++n ) std::clog << " ";
		std::clog << "] " << static_cast<int>(100*percent) << "%\r";
	}
}


// -------------------------------------------------------------------
// raytrace::print_location()
// �g���̈ʒu�A�g���x�N�g���A����єg�������܂�Ă���̌o�ߎ��Ԃ��o�͂���B
//
void raytrace::print_location( const raypath_element& ptr )
{
	const double Re = rtc::getCosmos().getPlanet().getRadius();
	const double t  = ptr.first;
	const rtc::vector r = ptr.second.first  / Re;// �ʒu�x�N�g��
	const rtc::vector k = ptr.second.second / Re;// �g���x�N�g��

	m_output
		<< t << " "
		<< r(0) << " " << r(1) << " " << r(2) << " "
		<< k(0) << " " << k(1) << " " << k(2) << "   ";
}

// -------------------------------------------------------------------
// raytrace::print_variation()
// �g���̈ʒu�A�g���x�N�g���A����єg�������܂�Ă���̌o�ߎ��Ԃ��o�͂���B
//
void raytrace::print_variation( const raypath_element& ptr )
{
	const double Re = rtc::getCosmos().getPlanet().getRadius();
	const double dt  = ptr.first;
	const rtc::vector dr = ptr.second.first  / Re;// �ʒu�x�N�g��
	const rtc::vector dk = ptr.second.second / Re;// �g���x�N�g��

	m_output
		<< dt << " "
		<< dr(0) << " " << dr(1) << " " << dr(2) << " "
		<< dk(0) << " " << dk(1) << " " << dk(2)
	<< "\n";
}

