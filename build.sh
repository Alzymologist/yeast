# !/bin/bash

if [ "$CF_PAGES_BRANCH" == "main" ]; then
    echo "Building for PRODUCTION 🔥"

    curl https://sh.rustup.rs -sSf | sh -s -- -y # Cargo installation
    source "$HOME/.cargo/env"
    # apt install graphviz # Graphviz installation
    
    # REPO_DIR="yeast_component"
    # REPO_URL="https://github.com/Alzymologist/yeast"
    # REPO_BRANCH="main"
  
    mkdir content
    cp -R content_fresh/* content 
    # mkdir $REPO_DIR
    # git clone --branch $REPO_BRANCH $REPO_URL $REPO_DIR
    # cd $REPO_DIR
    # cargo run -- --production
    # cp -R output/* ../static/yeast-component-output
    # cd .. 
    zola build

elif [ "$CF_PAGES_BRANCH" == "feature-main" ]; then
    echo "Building for STAGING 🔥"

    curl https://sh.rustup.rs -sSf | sh -s -- -y # Cargo installation
    source "$HOME/.cargo/env"
    apt install graphviz # Graphviz installation
    
    REPO_DIR="yeast_component"
    REPO_URL="https://github.com/Alzymologist/yeast"
    REPO_BRANCH="feature-2"
  
    mkdir content
    cp -R content_fresh/* content 
    mkdir $REPO_DIR
    git clone --branch $REPO_BRANCH $REPO_URL $REPO_DIR
    cd $REPO_DIR
    cargo run -- --staging
    cp -R output/* ../static/yeast-component-output
    cd .. 
    zola build

fi
