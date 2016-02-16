#ifndef BROWSEMODE_HPP
#define BROWSEMODE_HPP

#include "Mode.hpp"

#include <Character/Character.hpp>
#include <Character/Skin.hpp>
#include <Character/integrator.h>
#include <Character/System.h>

#include <vector>
#include <deque>
#include <utility>
#include <string>

using std::deque;
using std::vector;
using std::pair;
using std::string;

class BrowseMode : public Mode {
public:
	BrowseMode();
	virtual ~BrowseMode();

	virtual void update(float const elapsed_time);
	
	virtual void handle_event(SDL_Event const &event);

	virtual void draw();

	Vector3f camera;
	Vector3f target;
	bool track;

	Character::Pose current_pose;
	unsigned int current_motion;
	float time;
	float play_speed;
	bool dump_frames;
	int frame_number;
	
	unsigned char * buffer;

	Character::Skin skin;
	
	Vec3 obstacle_pos;
	Vec3 obstacle_vel;

	Integrator* integrator;
	System* cloth;
};

#endif //BROWSEMODE_HPP
