Yeast knowledge storage.

# OPTIONAL

Build script for VsCodium, which launches a component on Rust (processes data and populates website pages) and raises a local development server using Zola. If the build script is already running but the Zola server has not yet been launched, then relaunching the script will have no effect. If the Zola server is already running, then relaunching the script will stop the server, and another run will restart the build. The contents of the tasks.json file for VSCodium:   
```
{
  "version": "2.0.0",
  "tasks": [
    {
      "label": "local build",
      "type": "shell",
      "command": "./local_build.sh",
      "group": {
        "kind": "build",
      }
    },
    {
      "label": "stop zola",
      "type": "shell",
      "command": "pkill -f 'zola.*serve' || true",
      "presentation": {
        "reveal": "never"
      },
      "problemMatcher": []
    },
    {
      "label": "stop zola and local build",
      "dependsOn": ["stop zola", "local build"],
      "group": {
        "kind": "build",
        "isDefault": true
      },
      "problemMatcher": []
    }
  ]
}

```
