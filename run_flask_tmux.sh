#!/bin/bash

#####################################################################
#  FIRST TMUX SESSION — viraldb_flask
#####################################################################

SESSION_NAME_1="viraldb_flask"
CONDA_ENV_1="viraldb"
PROJECT_DIR_1="/mnt/c/Users/joeys/Documents/WorkingViralDB"
APP_CMD_1="python app.py"

tmux has-session -t $SESSION_NAME_1 2>/dev/null
if [ $? != 0 ]; then
    echo "🚀 Starting tmux session: $SESSION_NAME_1"
    tmux new-session -d -s $SESSION_NAME_1 \
        "source ~/miniconda3/etc/profile.d/conda.sh \
         && conda activate $CONDA_ENV_1 \
         && cd $PROJECT_DIR_1 \
         && $APP_CMD_1"
else
    echo "Session '$SESSION_NAME_1' already running."
fi



#####################################################################
#  SECOND TMUX SESSION — ligase_flask
#####################################################################

SESSION_NAME_2="ligase_flask"
CONDA_ENV_2="viraldb"
PROJECT_DIR_2="/mnt/c/Users/joeys/Documents/WorkingViralDB/Ligases/MODULE/e3-recruiter-mod"
APP_CMD_2="python Ligase_app.py"

tmux has-session -t $SESSION_NAME_2 2>/dev/null
if [ $? != 0 ]; then
    echo "🚀 Starting tmux session: $SESSION_NAME_2"
    tmux new-session -d -s $SESSION_NAME_2 \
        "source ~/miniconda3/etc/profile.d/conda.sh \
         && conda activate $CONDA_ENV_2 \
         && cd $PROJECT_DIR_2 \
         && $APP_CMD_2"
else
    echo "Session '$SESSION_NAME_2' already running."
fi



#####################################################################
#  THIRD TMUX SESSION — tts_flask  (NEW)
#####################################################################

SESSION_NAME_3="tts_flask"
CONDA_ENV_3="tts-flask"
PROJECT_DIR_3="/mnt/c/Users/joeys/Documents/TTS"
APP_CMD_3="python app.py"

tmux has-session -t $SESSION_NAME_3 2>/dev/null
if [ $? != 0 ]; then
    echo "🚀 Starting tmux session: $SESSION_NAME_3"
    tmux new-session -d -s $SESSION_NAME_3 \
        "source ~/miniconda3/etc/profile.d/conda.sh \
         && conda activate $CONDA_ENV_3 \
         && cd $PROJECT_DIR_3 \
         && $APP_CMD_3"
else
    echo "Session '$SESSION_NAME_3' already running."
fi



#####################################################################
# DONE
#####################################################################

echo "✔ All tmux Flask apps started (if not already running)."


# #!/bin/bash

# # Name of the tmux session
# SESSION_NAME="viraldb_flask"

# # Path to your Conda environment
# CONDA_ENV="viraldb"

# # Path to your project directory
# PROJECT_DIR="/mnt/c/Users/joeys/Documents/WorkingViralDB"

# # Command to run Flask app
# FLASK_CMD="python app.py"

# # Check if tmux session already exists
# tmux has-session -t $SESSION_NAME 2>/dev/null

# if [ $? != 0 ]; then
#     echo "Starting new tmux session: $SESSION_NAME"

#     # Start tmux and activate conda env, cd into dir, run Flask
#     tmux new-session -d -s $SESSION_NAME "source ~/miniconda3/etc/profile.d/conda.sh && conda activate $CONDA_ENV && cd $PROJECT_DIR && $FLASK_CMD"
# else
#     echo "Session $SESSION_NAME already running."
# fi

# # Attach to the session
# # tmux attach -t $SESSION_NAME



# #!/bin/bash

# #####################################################################
# #  FIRST TMUX SESSION — viraldb_flask
# #####################################################################

# SESSION_NAME_1="viraldb_flask"
# CONDA_ENV_1="viraldb"
# PROJECT_DIR_1="/mnt/c/Users/joeys/Documents/WorkingViralDB"
# APP_CMD_1="python app.py"

# tmux has-session -t $SESSION_NAME_1 2>/dev/null
# if [ $? != 0 ]; then
#     echo "🚀 Starting tmux session: $SESSION_NAME_1"
#     tmux new-session -d -s $SESSION_NAME_1 \
#         "source ~/miniconda3/etc/profile.d/conda.sh \
#          && conda activate $CONDA_ENV_1 \
#          && cd $PROJECT_DIR_1 \
#          && $APP_CMD_1"
# else
#     echo "Session '$SESSION_NAME_1' already running."
# fi



# #####################################################################
# #  SECOND TMUX SESSION — ligase_flask
# #####################################################################

# SESSION_NAME_2="ligase_flask"
# CONDA_ENV_2="viraldb"
# PROJECT_DIR_2="/mnt/c/Users/joeys/Documents/WorkingViralDB/Ligases/MODULE/e3-recruiter-mod"
# APP_CMD_2="python Ligase_app.py"

# tmux has-session -t $SESSION_NAME_2 2>/dev/null
# if [ $? != 0 ]; then
#     echo "🚀 Starting tmux session: $SESSION_NAME_2"
#     tmux new-session -d -s $SESSION_NAME_2 \
#         "source ~/miniconda3/etc/profile.d/conda.sh \
#          && conda activate $CONDA_ENV_2 \
#          && cd $PROJECT_DIR_2 \
#          && $APP_CMD_2"
# else
#     echo "Session '$SESSION_NAME_2' already running."
# fi


# #####################################################################
# # DONE
# #####################################################################

# echo "✔ All tmux Flask apps started (if not already running)."


